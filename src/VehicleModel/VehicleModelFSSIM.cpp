#include "VehicleModel/VehicleModelInterface.hpp"

#include "FSSIMVehicleModel.hpp"
#include <algorithm>
#include <cmath>

// Wrapper implementing IVehicleModel using the FSSIMVehicleModel via composition
class VehicleModelFSSIM : public IVehicleModel {
public:
    VehicleModelFSSIM() {
        this->position = Eigen::Vector3d::Zero();
        this->orientation = Eigen::Vector3d::Zero();
        this->velocity = Eigen::Vector3d::Zero();
        this->acceleration = Eigen::Vector3d::Zero();
        this->angularVelocity = Eigen::Vector3d::Zero();
        this->angularAcceleration = Eigen::Vector3d::Zero();

        this->torques = {0.0, 0.0, 0.0, 0.0};
        this->steeringAngles = {0.0, 0.0, 0.0, 0.0};
        this->wheelOrientations = {0.0, 0.0, 0.0, 0.0};
        this->wheelspeeds = {0.0, 0.0, 0.0, 0.0};

        // Reasonable defaults
        wheelRadius = 0.206;
        gearRatio = 12.0;
        innerSteeringRatio = 1.0;
        outerSteeringRatio = 1.0;
        nominalVoltageTS = 0.0;
    }

    // Map existing YAML from pacsim or leave defaults if not present
    bool readConfig(ConfigElement& config) override {
        try {
            auto cfg = config["simple_bicycle_model"];
            // Map basic geometry to FSSIM parameters if present
            double lr = fssim_.getParam().kinematic.l_R;
            double lf = fssim_.getParam().kinematic.l_F;
            double sf = fssim_.getParam().kinematic.b_F;
            double sr = fssim_.getParam().kinematic.b_R;

            cfg["kinematics"].getElement<double>(&lf, "lf");
            cfg["kinematics"].getElement<double>(&lr, "lr");
            cfg["kinematics"].getElement<double>(&sf, "sf");
            cfg["kinematics"].getElement<double>(&sr, "sr");

            // Update FSSIM parameters if provided
            auto& p = fssim_.params();
            p.kinematic.l = lf + lr;
            p.kinematic.b_F = sf;
            p.kinematic.b_R = sr;
            p.kinematic.w_front = (p.kinematic.l > 1e-6) ? (lr / p.kinematic.l) : p.kinematic.w_front;

            // Tire model mapping if present
            double Blat = p.tire.B, Clat = p.tire.C, Dlat = p.tire.D, Elat = p.tire.E;
            cfg["tire"].getElement<double>(&Blat, "Blat");
            cfg["tire"].getElement<double>(&Clat, "Clat");
            cfg["tire"].getElement<double>(&Dlat, "Dlat");
            cfg["tire"].getElement<double>(&Elat, "Elat");
            p.tire.tire_coefficient = 1.0;
            p.tire.B = Blat;
            p.tire.C = Clat;
            p.tire.D = Dlat;
            p.tire.E = Elat;

            // Aero map: convert to FSSIM constants
            double cla = 3.7, cda = 1.1, area = 1.1;
            cfg["aero"].getElement<double>(&cla, "cla");
            cfg["aero"].getElement<double>(&cda, "cda");
            cfg["aero"].getElement<double>(&area, "aeroArea");
            const double rho = 1.29;
            p.aero.c_down = 0.5 * rho * area * cla; // maps F_down = c_down * v^2
            p.aero.c_drag = 0.5 * rho * area * cda; // maps F_drag = c_drag * v^2

            // Mass/inertia
            cfg.getElement<double>(&p.inertia.m, "m");
            cfg.getElement<double>(&p.inertia.I_z, "Izz");

            // Drivetrain basics
            cfg.getElement<double>(&wheelRadius, "wheelRadius");
            cfg.getElement<double>(&gearRatio, "gearRatio");
            cfg.getElement<double>(&innerSteeringRatio, "innerSteeringRatio");
            cfg.getElement<double>(&outerSteeringRatio, "outerSteeringRatio");
            cfg.getElement<double>(&nominalVoltageTS, "nominalVoltageTS");

        // Derived params are recomputed internally when used
        return true;
        } catch (...) {
            // If config structure doesn't match, keep defaults
            return true;
        }
    }

    // Getters
    Eigen::Vector3d getPosition() override {
        const auto& s = fssim_.getState();
        this->position = Eigen::Vector3d(s.p.x(), s.p.y(), 0.0);
        return this->position;
    }
    Eigen::Vector3d getOrientation() override {
        const auto& s = fssim_.getState();
        this->orientation = Eigen::Vector3d(0.0, 0.0, s.yaw);
        return this->orientation;
    }
    Eigen::Vector3d getAngularVelocity() override { return this->angularVelocity; }
    Eigen::Vector3d getAngularAcceleration() override { return this->angularAcceleration; }
    Eigen::Vector3d getVelocity() override { return this->velocity; }
    Eigen::Vector3d getAcceleration() override { return this->acceleration; }
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
        Eigen::AngleAxisd yawAngle(this->orientation.z(), Eigen::Vector3d::UnitZ());
        auto R = yawAngle.matrix();
        Eigen::Vector3d FL = R * Eigen::Vector3d(lf, +sf * 0.5, 0.0) + this->position;
        Eigen::Vector3d FR = R * Eigen::Vector3d(lf, -sf * 0.5, 0.0) + this->position;
        Eigen::Vector3d RL = R * Eigen::Vector3d(-lr, +sr * 0.5, 0.0) + this->position;
        Eigen::Vector3d RR = R * Eigen::Vector3d(-lr, -sr * 0.5, 0.0) + this->position;
        return {FL, FR, RL, RR};
    }

    // Setters / inputs
    void setTorques(Wheels in) override { this->torques = in; }
    void setRpmSetpoints(Wheels in) override { this->rpmSetpoints = in; }
    void setMinTorques(Wheels in) override { this->minTorques = in; }
    void setMaxTorques(Wheels in) override { this->maxTorques = in; }
    void setSteeringSetpointFront(double in) override {
        // Map front steering command to symmetric front wheel angles and to FSSIM delta
        double avgRatio = 0.5 * (innerSteeringRatio + outerSteeringRatio);
        double delta = in * avgRatio;
        this->steeringAngles.FL = delta;
        this->steeringAngles.FR = delta;
        fssim_steer_ = delta;
    }
    void setSteeringSetpointRear(double) override { /* not supported */ }
    void setPowerGroundSetpoint(double in) override {
        // Treat as throttle [0..1]
        fssim_throttle_ = std::min(std::max(in, 0.0), 1.0);
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
        const auto prev = fssim_.getState();
        const double prev_r = prev.r;
        const Vec2 prev_v = prev.v;

        const auto& s = fssim_.step(fssim_throttle_, fssim_steer_, dt);

        // Update IVehicleModel state members
        this->position = Eigen::Vector3d(s.p.x(), s.p.y(), 0.0);
        this->orientation = Eigen::Vector3d(0.0, 0.0, s.yaw);
        this->velocity = Eigen::Vector3d(s.v.x(), s.v.y(), 0.0);
        this->angularVelocity = Eigen::Vector3d(0.0, 0.0, s.r);

        // Finite-difference accelerations in body frame
        const Vec2 a_body = (s.v - prev_v) / std::max(1e-6, dt);
        this->acceleration = Eigen::Vector3d(a_body.x(), a_body.y(), 0.0);
        this->angularAcceleration = Eigen::Vector3d(0.0, 0.0, (s.r - prev_r) / std::max(1e-6, dt));

        // Wheelspeeds: approximate from body x velocity and gear ratio
        double rpm = (s.v.x() / (wheelRadius * 2.0 * M_PI)) * 60.0 * gearRatio;
        this->wheelspeeds = {rpm, rpm, rpm, rpm};

        // Wheel orientations integrate trivially from wheelspeeds
        double dtheta = (rpm / (60.0 * gearRatio)) * dt * 2.0 * M_PI;
        this->wheelOrientations.FL = std::fmod(this->wheelOrientations.FL + dtheta, 2.0 * M_PI);
        this->wheelOrientations.FR = std::fmod(this->wheelOrientations.FR + dtheta, 2.0 * M_PI);
        this->wheelOrientations.RL = std::fmod(this->wheelOrientations.RL + dtheta, 2.0 * M_PI);
        this->wheelOrientations.RR = std::fmod(this->wheelOrientations.RR + dtheta, 2.0 * M_PI);
    }

private:
    // Underlying model
    FSSIMVehicleModel fssim_{};
    // Inputs
    double fssim_throttle_{0.0};
    double fssim_steer_{0.0};
    // Parameters needed for some interface methods
    double wheelRadius{0.206};
    double gearRatio{12.0};
    double innerSteeringRatio{1.0};
    double outerSteeringRatio{1.0};
    double nominalVoltageTS{0.0};

    Wheels minTorques{0,0,0,0};
    Wheels maxTorques{0,0,0,0};
    Wheels rpmSetpoints{0,0,0,0};
};
