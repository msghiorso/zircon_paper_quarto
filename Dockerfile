FROM registry.gitlab.com/enki-portal/thermoengine:master
COPY notebooks/1-Model-description.ipynb ${HOME}
COPY notebooks/2-Endmembers-MELTS.ipynb ${HOME}
COPY notebooks/3-Liquid-MELTS-codegen.ipynb ${HOME}
COPY notebooks/4-Liquid-MELTS-API-tests.ipynb ${HOME}
COPY notebooks/5-Liquid-MELTS-calib-1.ipynb ${HOME}
COPY notebooks/6-Liquid-MELTS-calib-2.ipynb ${HOME} 
COPY notebooks/7-MELTS-with-Zr.ipynb ${HOME}
COPY data/Boehnke_et_al.csv ${HOME}
COPY data/CO2_rMELTS_ZR_calib.c ${HOME}
COPY data/CO2_rMELTS_ZR_calib.h ${HOME}
COPY data/Gervasoni_2016.xlsx ${HOME}
COPY data/H2O_rMELTS_ZR_calib.c ${HOME}
COPY data/H2O_rMELTS_ZR_calib.h ${HOME}
COPY data/support_scripts.py ${HOME}
COPY data/Watson_and_Harrison.csv ${HOME}
COPY jupyter_notebook_config.py ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
