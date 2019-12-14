
# remote to local
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/03-covar_assoc/_rslurm_covar_corr /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/03-covar_assoc/
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/06-IPW_ML/_rslurm_vivid_spSL /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/06-IPW_ML/
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/06-IPW_ML/_rslurm_covar_IPTW_balance /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/06-IPW_ML/
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/_rslurm_prevmap_diagnostic_chains /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/07-spatial_prediction/
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_long_chains /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/07-spatial_prediction/
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/_rslurm_Prevmap_predictions /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/07-spatial_prediction/
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/09-Power/_rslurm_powercalc_pf /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/09-Power/
rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/09-Power/_rslurm_powercalc_pv /Users/nickbrazeau/Documents/GitHub/VivID_Epi/analyses/09-Power/

# local to remote
rsync -av /Users/nickbrazeau/Documents/GitHub/VivID_Epi/model_datamaps/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/model_datamaps/
rsync -av /Users/nickbrazeau/Documents/GitHub/VivID_Epi/data/derived_data nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/data/
rsync -av /Users/nickbrazeau/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/data/
rsync -avr /Users/nickbrazeau/Documents/GitHub/VivID_Epi/data/map_bases/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/data/