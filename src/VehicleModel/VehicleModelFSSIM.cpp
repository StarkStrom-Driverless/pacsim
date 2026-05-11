#include "VehicleModel/VehicleModelInterface.hpp"

#include "FSSIMVehicleModel.hpp"
#include <algorithm>
#include <cmath>

namespace {
constexpr double kCmdEps = 1e-6;
}

// Wrapper implementing IVehicleModel using the FSSIMVehicleModel via composition
class VehicleModelFSSIM : public IVehicleModel {
public:
    VehicleModelFSSIM() {
        // Reasonable defaults (kept minimal)
        this->torques = {0.0, 0.0, 0.0, 0.0, 0.0};
        this->steeringAngles = {0.0, 0.0, 0.0, 0.0, 0.0};
        this->wheelOrientations = {0.0, 0.0, 0.0, 0.0, 0.0};
        this->wheelspeeds = {0.0, 0.0, 0.0, 0.0, 0.0};

        wheelRadius = 0.206;
        gearRatio = 12.0;
        nominalVoltageTS = 0.0;
        base_g_ = fssim_.getParam().inertia.g;        
    }

    // Map existing YAML from pacsim or leave defaults if not present
    bool readConfig(ConfigElement& config) override {
        fssim_.loadCarConfig(config.getNode());
        return true;
    }

    // Getters
    Eigen::Vector3d getPosition() override {
        const auto& s = fssim_.getState();
        return Eigen::Vector3d(s.p.x(), s.p.y(), 0.0);
    }
    Eigen::Vector3d getOrientation() override {
        const auto& s = fssim_.getState();
        return Eigen::Vector3d(0.0, 0.0, s.yaw);
    }
    Eigen::Vector3d getAngularVelocity() override {
        const auto& s = fssim_.getState();
        return Eigen::Vector3d(0.0, 0.0, s.r);
    }
    Eigen::Vector3d getAngularAcceleration() override {
        if (last_dt_ <= 0.0) return Eigen::Vector3d::Zero();
        const auto& s = fssim_.getState();
        double rdot = (s.r - last_state_.r) / last_dt_;
        return Eigen::Vector3d(0.0, 0.0, rdot);
    }
    Eigen::Vector3d getVelocity() override {
        const auto& s = fssim_.getState();
        return Eigen::Vector3d(s.v.x(), s.v.y(), 0.0);
    }
    Eigen::Vector3d getAcceleration() override {
        const auto& s = fssim_.getState();
        return Eigen::Vector3d(s.a.x(), s.a.y(), 0.0);
    }
    Wheels getSteeringAngles() override { return this->steeringAngles; }
    Wheels getWheelspeeds() override { return this->wheelspeeds; }
    Wheels getTorques() override { return this->torques; }
    Wheels getWheelOrientations() override { return this->wheelOrientations; }
    double getSteeringWheelAngle() override { return this->steeringAngles.FL; }
    double getVoltageTS() override { return this->nominalVoltageTS; }
    double getCurrentTS() override { return 0.0; }

    std::array<Eigen::Vector3d, 4> getWheelPositions() override {
        const auto& p = fssim_.getParam();
        const double lf = p.kinematic.l_F;
        const double lr = p.kinematic.l_R;
        const double sf = p.kinematic.b_F;
        const double sr = p.kinematic.b_R;
        const auto& s = fssim_.getState();
        Eigen::AngleAxisd yawAngle(s.yaw, Eigen::Vector3d::UnitZ());
        auto R = yawAngle.matrix();
        Eigen::Vector3d base(s.p.x(), s.p.y(), 0.0);
        Eigen::Vector3d FL = R * Eigen::Vector3d(lf, +sf * 0.5, 0.0) + base;
        Eigen::Vector3d FR = R * Eigen::Vector3d(lf, -sf * 0.5, 0.0) + base;
        Eigen::Vector3d RL = R * Eigen::Vector3d(-lr, +sr * 0.5, 0.0) + base;
        Eigen::Vector3d RR = R * Eigen::Vector3d(-lr, -sr * 0.5, 0.0) + base;
        return {FL, FR, RL, RR};
    }

    // Setters / inputs
    void setTorques(Wheels in) override { this->torques = in; }
    void setRpmSetpoints(Wheels in) override {
        this->rpmSetpoints = in;
        this->hasRpmSetpoint_ = true;
        this->rpmCommandAgeS_ = 0.0;
    }
    void setMinTorques(Wheels in) override { this->minTorques = in; }

    void setMaxTorques(Wheels in) override { this->maxTorques = in; }
    void setSteeringSetpointFront(double in) override {
        // Single-track: front steer equals input
        double delta = in;
        this->steeringAngles.FL = delta;
        this->steeringAngles.FR = delta;
        fssim_steer_ = delta;
    }
    void setSteeringSetpointRear(double) override { /* not supported */ }
    void setPowerGroundSetpoint(double in) override {
        // Interpret command as longitudinal force in Newton from upstream controller.
        targetLongitudinalForceN_ = in;
        hasLongitudinalForceCommand_ = true;
        forceCommandAgeS_ = 0.0;
    }
    void setPosition(Eigen::Vector3d position) override {
        this->position = position;
        auto s = fssim_.getState();
        s.p.x() = position.x();
        s.p.y() = position.y();
        fssim_.setState(s);
    }
    void setOrientation(Eigen::Vector3d orientation) override {
        this->orientation = orientation;
        auto s = fssim_.getState();
        s.yaw = orientation.z();
        fssim_.setState(s);
    }

    void forwardIntegrate(double dt, Wheels /*frictionCoefficients*/) override {
        // Save previous for finite-difference accelerations
        const auto params = fssim_.getParam();
        const double total_max_torque = maxTorques.FL + maxTorques.FR + maxTorques.RL + maxTorques.RR;
        const bool use_force_command = hasLongitudinalForceCommand_;

        if (use_force_command) {
            // Map force command in N to drivetrain command input.
            fssim_dc_ = std::clamp(
                targetLongitudinalForceN_ / std::max(kCmdEps, params.driveTrain.cm1), -1.0, 1.0);

            // Optional safety: only apply rpm cap when the upstream actually provided rpm setpoints
            // (avoids accidental zeroing when rpmSetpoints are still at their default of 0).
            if (hasRpmSetpoint_) {
                double rpm_avg = 0.25 * (rpmSetpoints.FL + rpmSetpoints.FR + rpmSetpoints.RL + rpmSetpoints.RR);
                double omega = rpm_avg * 2.0 * M_PI / 60.0;
                double v_target = (omega * wheelRadius) / std::max(kCmdEps, gearRatio);
                if (v_target > kCmdEps) {
                    const double v_current = fssim_.getState().v.x();
                    if (v_current > v_target && fssim_dc_ > 0.0) {
                        fssim_dc_ = 0.0;
                    }
                }
            }
        } else if (std::abs(total_max_torque) > kCmdEps) {
            // Backward compatibility: legacy torque path.
            const double target_force = gearRatio * total_max_torque / std::max(kCmdEps, wheelRadius);
            fssim_dc_ = std::clamp(target_force / std::max(kCmdEps, params.driveTrain.cm1), -1.0, 1.0);
        } else {
            // Velocity path: use wheel speed setpoint as target body velocity.
            double rpm_avg = 0.25 * (rpmSetpoints.FL + rpmSetpoints.FR + rpmSetpoints.RL + rpmSetpoints.RR);
            double omega = rpm_avg * 2.0 * M_PI / 60.0;
            double v_target = (omega * wheelRadius) / std::max(kCmdEps, gearRatio);
            updateVelocityPI(v_target, dt);
        }

        // Save previous state for finite-difference getters
        last_state_ = fssim_.getState();

        const auto& s = fssim_.step(fssim_dc_, fssim_steer_, dt);

        last_dt_ = dt;

        // Wheelspeeds: approximate from body x velocity and gear ratio
        double rpm = (s.v.x() / (wheelRadius * 2.0 * M_PI)) * 60.0 * gearRatio;
        this->wheelspeeds = {rpm, rpm, rpm, rpm, 0.0};

        // Wheel orientations integrate trivially from wheelspeeds
        double dtheta = (rpm / (60.0 * gearRatio)) * dt * 2.0 * M_PI;
        this->wheelOrientations.FL = std::fmod(this->wheelOrientations.FL + dtheta, 2.0 * M_PI);
        this->wheelOrientations.FR = std::fmod(this->wheelOrientations.FR + dtheta, 2.0 * M_PI);
        this->wheelOrientations.RL = std::fmod(this->wheelOrientations.RL + dtheta, 2.0 * M_PI);
        this->wheelOrientations.RR = std::fmod(this->wheelOrientations.RR + dtheta, 2.0 * M_PI);

        // Keep the last force/rpm setpoints sticky so that PacSim ticks faster than the
        // upstream controller (e.g. 200 Hz sim vs 40 Hz MPC) do not fall back to the
        // velocity-PI brake-to-zero path between commands. Drop the command only after a
        // staleness timeout indicating the upstream publisher really stopped.
        constexpr double kCommandStaleTimeoutS = 0.2;
        if (hasLongitudinalForceCommand_) {
            forceCommandAgeS_ += dt;
            if (forceCommandAgeS_ > kCommandStaleTimeoutS) {
                hasLongitudinalForceCommand_ = false;
                targetLongitudinalForceN_ = 0.0;
                forceCommandAgeS_ = 0.0;
            }
        }
        if (hasRpmSetpoint_) {
            rpmCommandAgeS_ += dt;
            if (rpmCommandAgeS_ > kCommandStaleTimeoutS) {
                hasRpmSetpoint_ = false;
                rpmCommandAgeS_ = 0.0;
            }
        }
    }

private:
    // Simple PI controller to compute driver command based on target velocity
    void updateVelocityPI(double v_target, double dt) {
        auto params = fssim_.getParam();
        const auto& state = fssim_.getState();
        auto v_current = state.v.x();
        double error = v_target - v_current;

        vel_i_ += error * dt;
        auto p_anteil = params.velocity_controller.p * error;
        auto i_anteil =  params.velocity_controller.i * vel_i_;
        double target_dc = p_anteil + i_anteil;
        fssim_dc_ = std::max(-1.0, std::min(1.0, target_dc));
        //std::cout << "Velocity controller: target: " << v_target << ", error: " << error <<  "target_dc: " << target_dc   << ", p_anteil: " <<p_anteil << ", i_anteil: " << i_anteil << std::endl;
    }

    // Underlying model and state
    FSSIMVehicleModel fssim_{};
    double fssim_dc_{0.0};
    double fssim_steer_{0.0};
    double wheelRadius{0.206};
    double gearRatio{12.0};
    double nominalVoltageTS{0.0};
    double base_g_{9.81};
    double gravityFactor_{1.0};
    double targetLongitudinalForceN_{0.0};
    bool hasLongitudinalForceCommand_{false};
    bool hasRpmSetpoint_{false};
    double forceCommandAgeS_{0.0};
    double rpmCommandAgeS_{0.0};
    Wheels minTorques{ -0.0, -0.0, -0.0, -0.0, 0.0 };
    Wheels maxTorques{ 0.0, 0.0, 0.0, 0.0, 0.0 };
    Wheels rpmSetpoints{ 0.0, 0.0, 0.0, 0.0, 0.0 };
    FSSIMVehicleModel::State last_state_{};
    double last_dt_{0.0};
    
    
    double vel_i_{0.0};
};
