## Data analysis
### Data
Morphometric data are available at [DOI:10.18167/DVN1/PDPDS4](http://dx.doi.org/10.18167/DVN1/PDPDS4)
### Import data
```
data WORK.MORPHOMETRY    ;
	infile 'morphometry_apterous_mean_and_ratio.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
	informat specimen $14. ;
	informat MLL $1. ;
	informat species $8. ;
	format specimen $14. ;
	format MLL $1. ;
	format species $8. ;
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
