## Data analysis
### Data
Morphometric data are available at [DOI:10.18167/DVN1/PDPDS4](http://dx.doi.org/10.18167/DVN1/PDPDS4)
### Import data
```
data WORK.MORPHOMETRY    ;
  %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
  infile 'D:\Mes Données\Etudes\SCYLV\diversité vecteur\morphométrie & taxo\routines stats et data pour dépôt\morphometry_apterous_mean_and_ratio.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
  informat specimen $14. ;
  informat MLL $1. ;
  informat species $8. ;
  informat BL best32. ;
  informat cauda best32. ;
  informat caudaBW best32. ;
  informat urs best32. ;
  informat NsetaeCauda best32. ;
  informat htII best32. ;
  informat HindFemur best32. ;
  informat HindTibia best32. ;
  informat HindTibiaW best32. ;
  informat siph best32. ;
  informat siphDW best32. ;
  informat siphBW best32. ;
  informat AntI best32. ;
  informat AntII best32. ;
  informat AntIII_IV best32. ;
  informat AntV best32. ;
  informat VIb best32. ;
  informat pt best32. ;
  informat AntIIIBW best32. ;
  informat NsetaeAntIII_IV best32. ;
  informat NsetaeAntV best32. ;
  informat Ant best32. ;
  informat Ant_BL best32. ;
  informat urs_htII best32. ;
  informat pt_VIb best32. ;
  informat pt_cauda best32. ;
  informat pt_siph best32. ;
  informat HindTibia_pt best32. ;
  informat cauda_urs best32. ;
  informat urs_VIb best32. ;
  informat siph_BL best32. ;
  informat siph_siphBW best32. ;
  informat siph_cauda best32. ;
  informat cauda_caudaBW best32. ;
  format specimen $14. ;
  format MLL $1. ;
  format species $8. ;
  format BL best12. ;
  format cauda best12. ;
  format caudaBW best12. ;
  format urs best12. ;
  format NsetaeCauda best12. ;
  format htII best12. ;
  format HindFemur best12. ;
  format HindTibia best12. ;
  format HindTibiaW best12. ;
  format siph best12. ;
  format siphDW best12. ;
  format siphBW best12. ;
  format AntI best12. ;
  format AntII best12. ;
  format AntIII_IV best12. ;
  format AntV best12. ;
  format VIb best12. ;
  format pt best12. ;
  format AntIIIBW best12. ;
  format NsetaeAntIII_IV best12. ;
  format NsetaeAntV best12. ;
  format Ant best12. ;
  format Ant_BL best12. ;
  format urs_htII best12. ;
  format pt_VIb best12.; 
  format pt_cauda best12. ;
  format pt_siph best12. ;
  format HindTibia_pt best12. ;
  format cauda_urs best12. ;
  format urs_VIb best12. ;
  format siph_BL best12. ;
  format siph_siphBW best12. ;
  format siph_cauda best12. ;
  format cauda_caudaBW best12. ;
  input
    specimen $
    MLL $
    species $
    BL
    cauda
    caudaBW
    urs
    NsetaeCauda
    htII
    HindFemur
    HindTibia
    HindTibiaW
    siph
    siphDW
    siphBW
    AntI
    AntII
    AntIII_IV
    AntV
    VIb
    pt
    AntIIIBW
    NsetaeAntIII_IV
    NsetaeAntV
    Ant
    Ant_BL
    urs_htII
    pt_VIb
    pt_cauda
    pt_siph
    HindTibia_pt
    cauda_urs
    urs_VIb
    siph_BL
    siph_siphBW
    siph_cauda
    cauda_caudaBW
    ;
  if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
run;
```
### Data analysis
#### synthesis table (Supp. Table 1): computation of mean, min, max, by species
```
data morphometry_synth;
  set morphometry;
  * only specimens from MLL A, C, D, F and the Melanaphis sorghi paratype NHM-1915-81 are included in the computation;
  where MLL in("A", "C", "D", "F") or specimen="NHM-1915-81";
run;

proc sort data=morphometry_synth;
  by species;
run;
/* mean */
proc means data=morphometry_synth noprint;
  by species;
  output out=synth_mean mean=;
run;
/* minimum */
proc means data=morphometry_synth noprint;
  by species;
  output out=synth_min min=;
run;
/* maximum */
proc means data=morphometry_synth noprint;
  by species;
  output out=synth_max max=;
run;
/* observation number */
proc means data=morphometry_synth noprint;
  by species;
  output out=synth_nb n=;
run;
```
#### ANOVA
```
data result_pvalue;
  format dependent $15.;
run;
%macro compar(trait);
  ods output ModelANOVA = tata;
  proc glm data=morphometry_synth;
    class species;
    model &trait = species ;
  run; quit;
  data result_pvalue;
    set result_pvalue tata;
    if hypothesistype=3 then delete;
    keep dependent probf;
  run;
%mend;

%compar(BL);
%compar(htII);
%compar(HindFemur);
%compar(HindTibia);
%compar(HindTibiaW);
%compar(siph);
%compar(siphDW);
%compar(siphBW);
%compar(cauda);
%compar(caudaBW);
%compar(urs);
%compar(AntI);
%compar(AntII);
%compar(AntIII_IV);
%compar(AntV);
%compar(VIb);
%compar(pt);
%compar(ant);
%compar(AntIIIBW);
%compar(NsetaeCauda);
%compar(NsetaeAntIII_IV);
%compar(NsetaeAntV);
%compar(Ant_BL);
%compar(urs_htII);
%compar(pt_VIb);
%compar(pt_cauda);
%compar(pt_siph);
%compar(HindTibia_Pt);
%compar(cauda_urs);
%compar(urs_VIb);
%compar(siph_BL);
%compar(siph_siphBW);
%compar(siph_cauda);
%compar(cauda_caudaBW);

/* False Discovery Rate procedure to control the experiment-wise Type-1 error */
/* using Benjamini, Y. and Hochberg, Y. (1995) method                         */

data prob_fdr;
	set result_pvalue;
	where probf ne .;
	rename probf = Raw_P;
run;

proc multtest inpvalues=prob_fdr fdr out=corrected_p;
run;
```