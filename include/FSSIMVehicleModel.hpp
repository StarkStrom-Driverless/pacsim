/*
 * AMZ-Driverless
 * Copyright (c) 2018 Authors:
 *   - Juraj Kabzan <kabzanj@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// Minimal standalone FSSIM vehicle dynamics model copied from FSSIM
// Depends only on Eigen and yaml-cpp

#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <string>
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <iomanip>

// Convenience vector alias
using Vec2 = Eigen::Vector2d;

class FSSIMVehicleModel {
public:
    struct State {
        // Position in world frame and orientation
        Vec2   p{Vec2::Zero()};
        double yaw{0.0};

        // Body-frame velocities and yaw rate
        Vec2   v{Vec2::Zero()};
        double r{0.0};

        Vec2   a{Vec2::Zero()};

        State operator*(double dt) const {
            State s;
            s.p   = dt * p; s.yaw = dt * yaw;
            s.v   = dt * v; s.r = dt * r;
            s.a   = dt * a; 
            return s;
        }

        State operator+(const State& o) const {
            State s;
            s.p   = p + o.p; s.yaw = yaw + o.yaw;
            s.v   = v + o.v; s.r = r + o.r;
            s.a   = a + o.a;
            return s;
        }

        void validate() { }
    };

    struct Input {
        double dc{0.0};     // Throttle command in [-1, 1]
        double delta{0.0};  // Steering angle [rad]
    };

    struct Param {
        struct Inertia { double m{190.0}, g{9.81}, I_z{110.0}; } inertia;
        struct Kinematic {
            double l{1.53};     // wheelbase
            double b_F{1.22};   // front axle track width
            double b_R{1.22};   // rear axle track width
            double w_front{0.5};
            double l_F{0.0};    // derived: l*(1-w_front)
            double l_R{0.0};    // derived: l*w_front
        } kinematic;
        struct Tire { double tire_coefficient{1.0}, B{12.56}, C{-1.38}, D{1.60}, E{-0.58}; } tire;
        struct Aero { double c_down{1.22*2.6*0.6}, c_drag{0.7*1.0*1.0}; } aero;
        struct DriveTrain { double inertia{0.4}, r_dyn{0.231}, m_lon_add{0.0}, cm1{5000.0}, cr0{180.0}; int nm_wheels{4}; } driveTrain;
        struct VelocityController { double p{.1}, i{.1};} velocity_controller;
    };

public:
    FSSIMVehicleModel() { finalizeDerivedParams(); }



    bool loadCarConfig(const YAML::Node& cfg) {
        YAML::Node car = cfg["fssim_model"];
        if (!car) {
            std::cout << "WARNING: Could not load YAML" << std::endl;
            return false;
        }

        debug_ = car["debug_print"].as<bool>();
        // Inertia (only used fields)
        if (auto n = car["inertia"]) {
            param_.inertia.m        = n["m"].as<double>();
            param_.inertia.g        = n["g"].as<double>();
            param_.inertia.I_z      = n["I_z"].as<double>();
        }

        // Kinematics
        if (auto n = car["kinematics"]) {
            param_.kinematic.l       = n["l"].as<double>();
            param_.kinematic.b_F     = n["b_F"].as<double>();
            param_.kinematic.b_R     = n["b_R"].as<double>();
            param_.kinematic.w_front = n["w_front"].as<double>();
        }

        // Tire (Magic Formula parameters)
        if (auto n = car["tire"]) {
            param_.tire.tire_coefficient = n["tire_coefficient"].as<double>();
            param_.tire.B = n["B"].as<double>() / param_.tire.tire_coefficient;
            param_.tire.C = n["C"].as<double>();
            param_.tire.D = n["D"].as<double>() * param_.tire.tire_coefficient;
            param_.tire.E = n["E"].as<double>();
        }

        // Aero
        if (auto n = car["aero"]) {
            param_.aero.c_down = n["C_Down"]["a"].as<double>() * n["C_Down"]["b"].as<double>() * n["C_Down"]["c"].as<double>();
            param_.aero.c_drag = n["C_drag"]["a"].as<double>() * n["C_drag"]["b"].as<double>() * n["C_drag"]["c"].as<double>();
        }

        // Drivetrain
        if (auto n = car["drivetrain"]) {
            param_.driveTrain.cm1       = n["Cm1"].as<double>();
            param_.driveTrain.cr0       = n["Cr0"].as<double>();
            param_.driveTrain.inertia   = n["inertia"].as<double>();
            param_.driveTrain.r_dyn     = n["r_dyn"].as<double>();
            param_.driveTrain.nm_wheels = n["nm_wheels"].as<int>();
        }
        if (auto n = car["velocity_controller"]) {
            param_.velocity_controller.p       = n["p"].as<double>();
            param_.velocity_controller.i       = n["i"].as<double>();
        }

        finalizeDerivedParams();
        return true;
    }

    // Expose parameter reference for external configuration
    Param& params() { return param_; }
    const Param& getParam() const { return param_; }

    void setState(const State& s) { state_ = s; }
    const State& getState() const { return state_; }


    void setDebug(bool enabled) { debug_ = enabled; }

    // Integrates dynamics given dc (throttle), delta (steer), and time step dt [s]
    const State& step(double dc, double delta, double dt) {
        // First-order steering actuator model (simple low-pass)
        const double T_sample   = 0.1;
        const double T_steering = 0.145;
        input_.delta = (T_sample / (T_steering + T_sample)) * delta
                     + (T_steering / (T_steering + T_sample)) * old_delta_;
        old_delta_ = input_.delta;
        input_.dc    = dc;

        // Normal force incl. aero downforce
        const State  x_prev   = state_;
        const double Fz_total = getNormalForce(x_prev);

        // Lateral tire forces using bicycle model with split left/right
        double FyF_l, FyF_r, FyR_l, FyR_r;
        computeTireLateralForces(x_prev, input_, Fz_total, FyF_l, FyF_r, FyR_l, FyR_r);

        const double Fx   = getFx(x_prev, input_);

        const State x_dot = f(state_, input_, Fx, FyF_l + FyF_r, FyR_l + FyR_r, FyF_l, FyF_r);
        State x_next      = state_ + x_dot * dt;
        state_            = f_kin_correction(x_next, state_, input_, Fx, dt);
        state_.validate();

        // Instantaneous body-frame acceleration from forces (not via numerical derivative) (unused)
        /// TODO(Ivo) This is AI slop, is it correct ?
        {
            const double m_lon  = param_.inertia.m + param_.driveTrain.m_lon_add;
            const double FyFtot = FyF_l + FyF_r;
            const double FyRtot = FyR_l + FyR_r;
            const double v_x    = std::max(0.0, x_prev.v.x());
            const double ax     = (x_prev.r * x_prev.v.y()) + (Fx - std::sin(input_.delta) * (FyFtot)) / m_lon;
            const double ay     = ((std::cos(input_.delta) * FyFtot) + FyRtot) / param_.inertia.m - (x_prev.r * v_x);
            state_.a.x() = ax;
            state_.a.y() = ay;
        }

        if (debug_) {
            debugPrint(x_prev, Fx, Fz_total, FyF_l, FyF_r, FyR_l, FyR_r, state_);
        }
        return state_;
    }

private:
    void finalizeDerivedParams() {
        param_.kinematic.l_F   = param_.kinematic.l * (1.0 - param_.kinematic.w_front);
        param_.kinematic.l_R   = param_.kinematic.l * param_.kinematic.w_front;
        param_.driveTrain.m_lon_add = param_.driveTrain.nm_wheels * param_.driveTrain.inertia
                                      / (param_.driveTrain.r_dyn * param_.driveTrain.r_dyn);
    }

    // Magic Formula lateral tire force
    double pacejkaFy(double alpha, double Fz) const {
        const double B = param_.tire.B;
        const double C = param_.tire.C;
        const double D = param_.tire.D;
        const double E = param_.tire.E;
        const double mu_y = D * std::sin(C * std::atan(B * (1.0 - E) * alpha + E * std::atan(B * alpha)));
        return Fz * mu_y;
    }

    // Slip angles and lateral forces per wheel
    void computeTireLateralForces(const State& x,
                                  const Input& u,
                                  double Fz_total,
                                  double& FyF_l,
                                  double& FyF_r,
                                  double& FyR_l,
                                  double& FyR_r) const {
        const double v_x = std::max(1.0, x.v.x());

        // Front axle geometry
        const double l_F = param_.kinematic.l_F;
        const double b_F = param_.kinematic.b_F;
        // Rear axle geometry
        const double l_R = param_.kinematic.l_R;
        const double b_R = param_.kinematic.b_R;

        // Slip angles (left/right) for front and rear
        const double alphaF_l = std::atan((x.v.y() + (+1.0) * l_F * x.r) / (v_x - 0.5 * b_F * x.r)) - u.delta;
        const double alphaF_r = std::atan((x.v.y() + (+1.0) * l_F * x.r) / (v_x + 0.5 * b_F * x.r)) - u.delta;
        const double alphaR_l = std::atan((x.v.y() + (-1.0) * l_R * x.r) / (v_x - 0.5 * b_R * x.r));
        const double alphaR_r = std::atan((x.v.y() + (-1.0) * l_R * x.r) / (v_x + 0.5 * b_R * x.r));

        // Static load distribution to each wheel (includes 0.5 factor per wheel)
        const double Fz_front_wheel = 0.5 * (param_.kinematic.w_front) * Fz_total;
        const double Fz_rear_wheel  = 0.5 * (1.0 - param_.kinematic.w_front) * Fz_total;

        FyF_l = pacejkaFy(alphaF_l, Fz_front_wheel);
        FyF_r = pacejkaFy(alphaF_r, Fz_front_wheel);
        FyR_l = pacejkaFy(alphaR_l, Fz_rear_wheel);
        FyR_r = pacejkaFy(alphaR_r, Fz_rear_wheel);
    }

    // Continuous-time bicycle dynamics with split front lateral forces
    State f(const State& x,
            const Input& u,
            double Fx,
            double FyF_tot,
            double FyR_tot,
            double FyF_l,
            double FyF_r) const {
        const double v_x = std::max(0.0, x.v.x());
        const double m_lon = param_.inertia.m + param_.driveTrain.m_lon_add;

        State dx{};
        Eigen::Rotation2D<double> R(x.yaw);
        dx.p = R.toRotationMatrix() * x.v;
        dx.yaw = x.r;
        dx.v.x() = (x.r * x.v.y()) + (Fx - std::sin(u.delta) * (FyF_tot)) / m_lon;
        dx.v.y() = ((std::cos(u.delta) * FyF_tot) + FyR_tot) / param_.inertia.m - (x.r * v_x);
        dx.r   = ((std::cos(u.delta) * FyF_tot * param_.kinematic.l_F
                  + std::sin(u.delta) * (FyF_l - FyF_r) * 0.5 * param_.kinematic.b_F)
                 - (FyR_tot * param_.kinematic.l_R)) / param_.inertia.I_z;
        dx.a.setZero();
        return dx;
    }

    // Low-speed kinematic correction and blending
    State f_kin_correction(const State& x_in,
                           const State& x_prev,
                           const Input& u,
                           double Fx,
                           double dt) const {
        State x = x_in;
        const double v_x_dot = Fx / (param_.inertia.m + param_.driveTrain.m_lon_add);
        const double v       = x_prev.v.norm();
        const double v_blend = 0.5 * (v - 1.5);
        const double blend   = std::max(0.0, std::min(1.0, v_blend));

        x.v.x() = blend * x.v.x() + (1.0 - blend) * (x_prev.v.x() + dt * v_x_dot);

        const double v_y_kin = std::tan(u.delta) * x.v.x() * param_.kinematic.l_R / param_.kinematic.l;
        const double r_kin   = std::tan(u.delta) * x.v.x() / param_.kinematic.l;

        x.v.y() = blend * x.v.y() + (1.0 - blend) * v_y_kin;
        x.r   = blend * x.r   + (1.0 - blend) * r_kin;
        return x;
    }

    // Aerodynamics + gravity normal force approximation
    double getNormalForce(const State& x) const {
        return param_.inertia.g * param_.inertia.m + param_.aero.c_down * x.v.x() * x.v.x();
    }

    // Longitudinal force from throttle minus drag and rolling resistance
    double getFx(const State& x, const Input& u) const {
        double F_drag = param_.aero.c_drag * x.v.x() * x.v.x();
        if(std::abs(F_drag) < 0.1) 
            F_drag = 0;
        
        auto rollwiderstand = param_.driveTrain.cr0;
        if(x.v.x() < 0.05)  
            rollwiderstand = 0;
        return dc * param_.driveTrain.cm1 - F_drag - rollwiderstand;
    }

private:
    Param param_{};
    State state_{};
    Input input_{};
    double old_delta_{0.0};
    bool debug_{false};

    void debugPrint(const State& x_prev,
                    double Fx,
                    double Fz,
                    double FyF_l,
                    double FyF_r,
                    double FyR_l,
                    double FyR_r,
                    const State& x_next) const {

        std::cout << std::setprecision(4);
        std::cout << "[FSSIM] step()" << std::endl;
        std::cout << " input: " << std::endl;
        std::cout << "  dc:    " << input_.dc << std::endl;
        std::cout << "  delta: " << input_.delta << std::endl;
        std::cout << " state (prev):" << std::endl;
        std::cout << "  p.x:  " << x_prev.p.x() << std::endl;
        std::cout << "  p.y:  " << x_prev.p.y() << std::endl;
        std::cout << "  yaw:  " << x_prev.yaw << std::endl;
        std::cout << "  v.x:  " << x_prev.v.x() << std::endl;
        std::cout << "  v.y:  " << x_prev.v.y() << std::endl;
        std::cout << "  r:    " << x_prev.r << std::endl;
        std::cout << " forces:" << std::endl;
        std::cout << "  Fz:   " << Fz << std::endl;
        std::cout << "  Fx:   " << Fx << std::endl;
        std::cout << "  FyF_l:" << FyF_l << std::endl;
        std::cout << "  FyF_r:" << FyF_r << std::endl;
        std::cout << "  FyR_l:" << FyR_l << std::endl;
        std::cout << "  FyR_r:" << FyR_r << std::endl;
        std::cout << " state (next):" << std::endl;
        std::cout << "  p.x:  " << x_next.p.x() << std::endl;
        std::cout << "  p.y:  " << x_next.p.y() << std::endl;
        std::cout << "  yaw:  " << x_next.yaw << std::endl;
        std::cout << "  v.x:  " << x_next.v.x() << std::endl;
        std::cout << "  v.y:  " << x_next.v.y() << std::endl;
        std::cout << "  r:    " << x_next.r << std::endl;
        std::cout << "  a.x:  " << x_next.a.x() << std::endl;
        std::cout << "  a.y:  " << x_next.a.y() << std::endl;
    }
};
