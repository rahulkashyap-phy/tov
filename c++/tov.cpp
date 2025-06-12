#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <chrono>
#include <numeric>

// Physical constants (CGS units)
namespace Constants {
    const double CGS_G = 6.67384e-8;      // cm^3 g^-1 s^-2
    const double CGS_C = 29979245800.0;   // cm s^-1
    const double CGS_MSUN = 1.98855e33;   // g
    const double CGS_AMU = 1.660538921e-24; // g
    const double fm = 1.0e-13;            // femtometer in cm
    const double dens_conversion = CGS_AMU / (fm * fm * fm);
    const double edens_conversion = CGS_C * CGS_C;
}

// Structure to hold EOS data
struct EOSData {
    std::vector<double> lg_rho;
    std::vector<double> lg_edens;
    std::vector<double> lg_pres;
    std::vector<double> gamma_arr;
};

// Structure to hold TOV solution state
struct TOVState {
    double P;        // Pressure
    double m;        // Mass
    double m_baryon; // Baryon mass
    double yp;       // Y potential
};

// Linear interpolation function
double interpolate(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    if (xi <= x[0]) return y[0];
    if (xi >= x.back()) return y.back();
    
    auto it = std::lower_bound(x.begin(), x.end(), xi);
    size_t i = std::distance(x.begin(), it) - 1;
    
    if (i >= x.size() - 1) return y.back();
    
    double t = (xi - x[i]) / (x[i+1] - x[i]);
    return y[i] * (1.0 - t) + y[i+1] * t;
}

// Calculate gradient using central differences
std::vector<double> gradient(const std::vector<double>& y, const std::vector<double>& x) {
    std::vector<double> grad(y.size());
    
    // Forward difference for first point
    grad[0] = (y[1] - y[0]) / (x[1] - x[0]);
    
    // Central difference for interior points
    for (size_t i = 1; i < y.size() - 1; ++i) {
        grad[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1]);
    }
    
    // Backward difference for last point
    size_t n = y.size() - 1;
    grad[n] = (y[n] - y[n-1]) / (x[n] - x[n-1]);
    
    return grad;
}

class TOVSolver {
private:
    EOSData eos_data;
    std::string eos_file;
    
public:
    TOVSolver(const std::string& eos_filename) : eos_file(eos_filename) {
        loadEOSData();
    }
    
    void loadEOSData() {
        std::ifstream file(eos_file);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open EOS file: " + eos_file);
        }
        
        std::string line;
        int skip_lines = 9;
        for (int i = 0; i < skip_lines; ++i) {
            std::getline(file, line);
        }
        
        std::vector<double> rho, edens, pres;
        double index, val1, val2, val3;
        
        while (file >> index >> val1 >> val2 >> val3) {
            rho.push_back(val1);
            edens.push_back(val2);
            pres.push_back(val3);
        }
        file.close();
        
        // Convert to log scale
        eos_data.lg_rho.resize(rho.size());
        eos_data.lg_edens.resize(edens.size());
        eos_data.lg_pres.resize(pres.size());
        
        for (size_t i = 0; i < rho.size(); ++i) {
            eos_data.lg_rho[i] = std::log10(rho[i]);
            eos_data.lg_edens[i] = std::log10(edens[i]);
            eos_data.lg_pres[i] = std::log10(pres[i]);
        }
        
        // Sort by density
        std::vector<size_t> indices(rho.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                  [&](size_t i, size_t j) { return eos_data.lg_rho[i] < eos_data.lg_rho[j]; });
        
        std::vector<double> sorted_lg_rho(rho.size());
        std::vector<double> sorted_lg_edens(edens.size());
        std::vector<double> sorted_lg_pres(pres.size());
        
        for (size_t i = 0; i < indices.size(); ++i) {
            sorted_lg_rho[i] = eos_data.lg_rho[indices[i]];
            sorted_lg_edens[i] = eos_data.lg_edens[indices[i]];
            sorted_lg_pres[i] = eos_data.lg_pres[indices[i]];
        }
        
        eos_data.lg_rho = sorted_lg_rho;
        eos_data.lg_edens = sorted_lg_edens;
        eos_data.lg_pres = sorted_lg_pres;
        
        // Calculate gamma array
        eos_data.gamma_arr = gradient(eos_data.lg_pres, eos_data.lg_edens);
    }
    
    // Get EOS quantities from pressure
    void eosFromPres(double pres_in, double& rho_out, double& edens_out, 
                     double& gamma_out, double& cs_out) {
        double lg_pres_in = std::log10(pres_in);
        
        double lg_rho_interp = interpolate(eos_data.lg_pres, eos_data.lg_rho, lg_pres_in);
        double lg_edens_interp = interpolate(eos_data.lg_pres, eos_data.lg_edens, lg_pres_in);
        double gamma_interp = interpolate(eos_data.lg_pres, eos_data.gamma_arr, lg_pres_in);
        
        rho_out = std::pow(10.0, lg_rho_interp) * Constants::dens_conversion;
        edens_out = std::pow(10.0, lg_edens_interp);
        gamma_out = gamma_interp;
        cs_out = std::sqrt((pres_in / edens_out) * gamma_interp);
    }
    
    // Get EOS quantities from density
    void eosFromDens(double rho_in, double& pres_out, double& edens_out) {
        double lg_rho_in = std::log10(rho_in / Constants::dens_conversion);
        
        double lg_pres_interp = interpolate(eos_data.lg_rho, eos_data.lg_pres, lg_rho_in);
        double lg_edens_interp = interpolate(eos_data.lg_rho, eos_data.lg_edens, lg_rho_in);
        
        pres_out = std::pow(10.0, lg_pres_interp);
        edens_out = std::pow(10.0, lg_edens_interp);
    }
    
    // TOV equations
    std::vector<double> tovEquations(double r, const std::vector<double>& y) {
        double P = y[0];
        double m = y[1];
        double m_baryon = y[2];
        double yp = y[3];
        
        double rho, eden, gamma, cs;
        eosFromPres(P, rho, eden, gamma, cs);
        
        double G = Constants::CGS_G;
        double c = Constants::CGS_C;
        
        double dPdr = -G * (eden + P / (c * c)) * (m + 4.0 * M_PI * r * r * r * P / (c * c));
        dPdr = dPdr / (r * (r - 2.0 * G * m / (c * c)));
        
        double dmdr = 4.0 * M_PI * r * r * eden;
        
        double dm_baryondr = dmdr / std::sqrt(1.0 - 2.0 * G * m / (r * c * c));
        
        rho = eden * Constants::CGS_C * Constants::CGS_C;
        double dypdr = -yp * yp / r - 
                       (r + (G / (c * c * c * c)) * 4.0 * M_PI * r * r * r * (P - rho)) * yp / 
                       (r * (r - 2.0 * G * m / (c * c))) +
                       (G * G / (c * c * c * c)) * 
                       (4.0 * (m + 4.0 * M_PI * r * r * r * P / (c * c)) * 
                        (m + 4.0 * M_PI * r * r * r * P / (c * c))) / 
                       (r * (r - 2.0 * G * m / (c * c)) * (r - 2.0 * G * m / (c * c))) +
                       6.0 / (r - 2.0 * Constants::CGS_G * m / Constants::CGS_C / Constants::CGS_C) -
                       4.0 * M_PI * (r * r) * (5.0 * rho + 9.0 * P + (rho + P) * (rho + P) / (P * gamma)) *
                       G / (c * c * c * c * (r - 2.0 * G * m / (c * c)));
        
        return {dPdr, dmdr, dm_baryondr, dypdr};
    }
    
    // Simple RK4 integrator
    std::vector<std::vector<double>> integrate(double r_start, double r_end, 
                                               const std::vector<double>& y0, int num_points) {
        std::vector<std::vector<double>> result(4);
        std::vector<double> r_array(num_points);
        
        double h = (r_end - r_start) / (num_points - 1);
        std::vector<double> y = y0;
        
        for (int i = 0; i < num_points; ++i) {
            double r = r_start + i * h;
            r_array[i] = r;
            
            for (int j = 0; j < 4; ++j) {
                result[j].push_back(y[j]);
            }
            
            if (i < num_points - 1) {
                // RK4 step
                std::vector<double> k1 = tovEquations(r, y);
                
                std::vector<double> y_temp(4);
                for (int j = 0; j < 4; ++j) {
                    y_temp[j] = y[j] + 0.5 * h * k1[j];
                }
                std::vector<double> k2 = tovEquations(r + 0.5 * h, y_temp);
                
                for (int j = 0; j < 4; ++j) {
                    y_temp[j] = y[j] + 0.5 * h * k2[j];
                }
                std::vector<double> k3 = tovEquations(r + 0.5 * h, y_temp);
                
                for (int j = 0; j < 4; ++j) {
                    y_temp[j] = y[j] + h * k3[j];
                }
                std::vector<double> k4 = tovEquations(r + h, y_temp);
                
                for (int j = 0; j < 4; ++j) {
                    y[j] += (h / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
                }
            }
        }
        
        return result;
    }
    
    // Find surface using bisection method
    double findSurface(double pmin, double rhoc, double rad_high) {
        auto surface_func = [this, pmin, rhoc](double r_max) -> double {
            double pres_c, edens_c;
            eosFromDens(rhoc, pres_c, edens_c);
            
            double r_min = 1.0e-3;
            double r3 = r_min * r_min * r_min;
            double m = (4.0 / 3.0) * M_PI * r3 * edens_c;
            double m_baryon = m / std::sqrt(1.0 - 2.0 * Constants::CGS_G * m / (r_min * Constants::CGS_C * Constants::CGS_C));
            double yp = 2.0;
            
            std::vector<double> y0 = {pres_c, m, m_baryon, yp};
            auto result = integrate(r_min, r_max, y0, 2000);
            
            return result[0].back() - pmin;
        };
        
        // Bisection method
        double a = 6.0e5, b = rad_high;
        double tolerance = 1.0e-6;
        
        while (std::abs(b - a) > tolerance) {
            double c = (a + b) / 2.0;
            if (surface_func(c) * surface_func(a) < 0) {
                b = c;
            } else {
                a = c;
            }
        }
        
        return (a + b) / 2.0;
    }
    
    // Calculate tidal deformability
    double calcTidalDeformability(double C, double Y) {
        double zeta = 4.0 * C * C * C * (13.0 - 11.0 * Y + C * (3.0 * Y - 2.0) + 
                     2.0 * (C * C) * (1.0 + Y)) +
                     3.0 * ((1.0 - 2.0 * C) * (1.0 - 2.0 * C)) * 
                     (2.0 - Y + 2.0 * C * (Y - 1.0)) * std::log(1.0 - 2.0 * C) +
                     2.0 * C * (6.0 - 3.0 * Y + 3.0 * C * (5.0 * Y - 8.0));
        
        double Lambda_dimensionless = (16.0 / (15.0 * zeta)) * 
                                      ((1.0 - 2.0 * C) * (1.0 - 2.0 * C)) * 
                                      (2.0 + 2.0 * C * (Y - 1.0) - Y);
        
        return Lambda_dimensionless;
    }
    
    // Main TOV solver
    void solve() {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::cout << "EOS file: " << eos_file << std::endl;
        
        // Extract EOS key from filename
        std::string eos_key = eos_file;
        size_t pos1 = eos_key.find("./eos_tables/");
        if (pos1 != std::string::npos) {
            eos_key = eos_key.substr(pos1 + 13);
        }
        size_t pos2 = eos_key.find(".lorene");
        if (pos2 != std::string::npos) {
            eos_key = eos_key.substr(0, pos2);
        }
        
        // Open output file
        std::ofstream outfile("tov_" + eos_key + ".txt");
        outfile << "#for eos = " << eos_file << "\n";
        outfile << "#Grav_Mass (solar mass) Radius (km) Lambda (dimensionless) rhoc (gm/cm^3) Compactness (dimensionless) Baryon_Mass (solar mass)\n";
        
        double pmin = 1.0e-10;
        int len_seq = 5;
        
        // Generate rhoc array
        std::vector<double> rhoc_arr(len_seq);
        double log_start = std::log10(6.0e14);
        double log_end = std::log10(5.0e15);
        for (int i = 0; i < len_seq; ++i) {
            double log_val = log_start + i * (log_end - log_start) / (len_seq - 1);
            rhoc_arr[i] = std::pow(10.0, log_val);
        }
        
        for (double rhoc : rhoc_arr) {
            std::cout << "rhoc = " << rhoc << std::endl;
            
            // Debug: Get initial EOS values
            double pres_c, edens_c;
            eosFromDens(rhoc, pres_c, edens_c);
            std::cout << "Initial: P_c = " << pres_c << ", edens_c = " << edens_c << std::endl;

            auto step_start = std::chrono::high_resolution_clock::now();
            
            // Find stellar radius
            double rstar = findSurface(pmin, rhoc, 3.0e6);
            
            auto step_end = std::chrono::high_resolution_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start);
            std::cout << "Time elapsed in root finding: " << step_duration.count() << " ms" << std::endl;
            
            // Solve TOV equations
            double r_min = 1.0e-3;
            double r3 = r_min * r_min * r_min;
            double m = (4.0 / 3.0) * M_PI * r3 * edens_c;
            double m_baryon = m / std::sqrt(1.0 - 2.0 * Constants::CGS_G * m / (r_min * Constants::CGS_C * Constants::CGS_C));
            double yp = 2.0;
            
            std::cout << "Initial conditions: m = " << m << ", m_baryon = " << m_baryon << ", yp = " << yp << std::endl;
            
            std::vector<double> y0 = {pres_c, m, m_baryon, yp};
            auto result = integrate(r_min, rstar, y0, 2000);
            
            double final_P = result[0].back();
            double final_m = result[1].back();
            double final_m_baryon = result[2].back();
            double final_yp = result[3].back();
            
            std::cout << std::fixed << std::setprecision(2)
                      << "FINAL: rstar: " << rstar / 1.0e5 << ", grav. mass: " << final_m / Constants::CGS_MSUN
                      << ", bary. mass: " << final_m_baryon / Constants::CGS_MSUN << ", yp: " << final_yp << std::endl;
            
            double C = (Constants::CGS_G / (Constants::CGS_C * Constants::CGS_C)) * final_m / rstar;
            double lambda_dimensionless = calcTidalDeformability(C, final_yp);
            
            outfile << std::fixed << std::setprecision(4)
                    << final_m / Constants::CGS_MSUN << "  "
                    << rstar / 1.0e5 << "  "
                    << lambda_dimensionless << "  "
                    << std::scientific << rhoc << "  "
                    << std::fixed << C << "  "
                    << final_m_baryon / Constants::CGS_MSUN << "\n";
            
            std::cout << std::fixed << std::setprecision(4)
                      << "FINAL: rstar: " << rstar / 1.0e5 << ", grav. mass: " << final_m / Constants::CGS_MSUN
                      << ", bary. mass: " << final_m_baryon / Constants::CGS_MSUN << ", yp: " << final_yp
                      << ", lambda: " << lambda_dimensionless << ", rhoc: " << std::scientific << rhoc
                      << ", compactness: " << std::fixed << C << std::endl;
        }
        
        outfile.close();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        std::cout << "Total time elapsed: " << total_duration.count() << " seconds" << std::endl;
    }
};

int main() {
    try {
        std::string eos_file = "../eos_tables/eosSLy.lorene";
        TOVSolver solver(eos_file);
        solver.solve();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

