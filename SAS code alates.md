## Data analysis
### Data
Morphometric data are available at [DOI:10.18167/DVN1/PDPDS4](http://dx.doi.org/10.18167/DVN1/PDPDS4)
### Import data
```
data WORK.morphometry    ;
  infile 'morphometry_alates_mean_and_ratio.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
  informat specimen $15. ;
  informat species $8. ;
  format specimen $15. ;
  format species $8. ;
  input
    specimen $
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
    NRhinAntIII
    NRhinAntIV
    NRhinAntV
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
run;
```
### Data analysis
#### synthesis table (Supp. Table 2): computation of mean, min, max, by species
```
data morphometry_synth;
  set morphometry;
  * discard the alate specimen from Theobald's type series which is low quality;
  where specimen ne "BM-1930-204-1";
run;

proc sort data=morphometry_synth ;
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

