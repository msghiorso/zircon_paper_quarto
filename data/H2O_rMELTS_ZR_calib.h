
const char *H2O_rMELTS_ZR_calib_identifier(void);
const char *H2O_rMELTS_ZR_calib_name(void);
const char *H2O_rMELTS_ZR_calib_formula(void);
const double H2O_rMELTS_ZR_calib_mw(void);
const double *H2O_rMELTS_ZR_calib_elements(void);

double H2O_rMELTS_ZR_calib_g(double T, double P);
double H2O_rMELTS_ZR_calib_dgdt(double T, double P);
double H2O_rMELTS_ZR_calib_dgdp(double T, double P);
double H2O_rMELTS_ZR_calib_d2gdt2(double T, double P);
double H2O_rMELTS_ZR_calib_d2gdtdp(double T, double P);
double H2O_rMELTS_ZR_calib_d2gdp2(double T, double P);
double H2O_rMELTS_ZR_calib_d3gdt3(double T, double P);
double H2O_rMELTS_ZR_calib_d3gdt2dp(double T, double P);
double H2O_rMELTS_ZR_calib_d3gdtdp2(double T, double P);
double H2O_rMELTS_ZR_calib_d3gdp3(double T, double P);

double H2O_rMELTS_ZR_calib_s(double T, double P);
double H2O_rMELTS_ZR_calib_v(double T, double P);
double H2O_rMELTS_ZR_calib_cv(double T, double P);
double H2O_rMELTS_ZR_calib_cp(double T, double P);
double H2O_rMELTS_ZR_calib_dcpdt(double T, double P);
double H2O_rMELTS_ZR_calib_alpha(double T, double P);
double H2O_rMELTS_ZR_calib_beta(double T, double P);
double H2O_rMELTS_ZR_calib_K(double T, double P);
double H2O_rMELTS_ZR_calib_Kp(double T, double P);

int H2O_rMELTS_ZR_get_param_number(void);
const char **H2O_rMELTS_ZR_get_param_names(void);
const char **H2O_rMELTS_ZR_get_param_units(void);
void H2O_rMELTS_ZR_get_param_values(double **values);
int H2O_rMELTS_ZR_set_param_values(double *values);
double H2O_rMELTS_ZR_get_param_value(int index);
int H2O_rMELTS_ZR_set_param_value(int index, double value);

double H2O_rMELTS_ZR_dparam_g(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_dgdt(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_dgdp(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_d2gdt2(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_d2gdtdp(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_d2gdp2(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_d3gdt3(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_d3gdt2dp(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_d3gdtdp2(double T, double P, int index);
double H2O_rMELTS_ZR_dparam_d3gdp3(double T, double P, int index);
