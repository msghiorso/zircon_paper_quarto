
const char *CO2_rMELTS_ZR_calib_identifier(void);
const char *CO2_rMELTS_ZR_calib_name(void);
const char *CO2_rMELTS_ZR_calib_formula(void);
const double CO2_rMELTS_ZR_calib_mw(void);
const double *CO2_rMELTS_ZR_calib_elements(void);

double CO2_rMELTS_ZR_calib_g(double T, double P);
double CO2_rMELTS_ZR_calib_dgdt(double T, double P);
double CO2_rMELTS_ZR_calib_dgdp(double T, double P);
double CO2_rMELTS_ZR_calib_d2gdt2(double T, double P);
double CO2_rMELTS_ZR_calib_d2gdtdp(double T, double P);
double CO2_rMELTS_ZR_calib_d2gdp2(double T, double P);
double CO2_rMELTS_ZR_calib_d3gdt3(double T, double P);
double CO2_rMELTS_ZR_calib_d3gdt2dp(double T, double P);
double CO2_rMELTS_ZR_calib_d3gdtdp2(double T, double P);
double CO2_rMELTS_ZR_calib_d3gdp3(double T, double P);

double CO2_rMELTS_ZR_calib_s(double T, double P);
double CO2_rMELTS_ZR_calib_v(double T, double P);
double CO2_rMELTS_ZR_calib_cv(double T, double P);
double CO2_rMELTS_ZR_calib_cp(double T, double P);
double CO2_rMELTS_ZR_calib_dcpdt(double T, double P);
double CO2_rMELTS_ZR_calib_alpha(double T, double P);
double CO2_rMELTS_ZR_calib_beta(double T, double P);
double CO2_rMELTS_ZR_calib_K(double T, double P);
double CO2_rMELTS_ZR_calib_Kp(double T, double P);

int CO2_rMELTS_ZR_get_param_number(void);
const char **CO2_rMELTS_ZR_get_param_names(void);
const char **CO2_rMELTS_ZR_get_param_units(void);
void CO2_rMELTS_ZR_get_param_values(double **values);
int CO2_rMELTS_ZR_set_param_values(double *values);
double CO2_rMELTS_ZR_get_param_value(int index);
int CO2_rMELTS_ZR_set_param_value(int index, double value);

double CO2_rMELTS_ZR_dparam_g(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_dgdt(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_dgdp(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_d2gdt2(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_d2gdtdp(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_d2gdp2(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_d3gdt3(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_d3gdt2dp(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_d3gdtdp2(double T, double P, int index);
double CO2_rMELTS_ZR_dparam_d3gdp3(double T, double P, int index);
