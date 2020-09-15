## Data analysis
### Data
Morphometric data are available at [DOI:10.18167/DVN1/PDPDS4](http://dx.doi.org/10.18167/DVN1/PDPDS4)
### Import data
```
data WORK.morphometry;
	infile 'morphometry_apterous_raw.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2;
	informat specimen $14.;
	informat MLL $1.;
	informat species $8.;
	format specimen $14.;
	format MLL $1.;
	input
	specimen $
	MLL $
	species $
	BL
	htII_r
	htII_l
	HindFemur_r
	HindFemur_l
	HindTibia_r
	HindTibia_l
	HindTibiaW_r
	HindTibiaW_l
	siph_r
	siphDW_r
	siphBW_r
	siph_l
	siphDW_l
	siphBW_l
	cauda
	caudaBW
	urs
	AntI_l
	AntII_l
	AntIII_l
	AntIIIBW_l
	AntIV_l
	AntV_l
	VIb_l
	pt_l
	Ant_l
	AntI_r
	AntII_r
	AntIII_r
	AntIIIBW_r
	AntIV_r
	AntV_r
	VIb_r
	pt_r
	Ant_r
	NsetaeCauda
	NsetaeAntIII_IV_r
	NsetaeAntIII_IV_l $
	NsetaeAntV_r
	NsetaeAntV_l $
	;
run;
```
### Computation of means (right and left sides), and ratio
```
data morpho2;
	set morphometry;
	antIII_IV_r=sum(of antiii_r, antiv_r);
	antIII_IV_l=sum(of antiii_l, antiv_l);
	htII=mean(of htII_r, htII_l);
	HindFemur=mean(of HindFemur_r, HindFemur_l);
	HindTibia=mean(of HindTibia_r, HindTibia_l);
	HindTibiaW=mean(of HindTibiaW_r, HindTibiaW_l);
	siph=mean(of siph_r, siph_l);
	siphDW=mean(of siphDW_r, siphDW_l);
	siphBW=mean(of siphBW_r, siphBW_l);
	AntI=mean(of AntI_r, AntI_l);
	AntII=mean(of AntII_r, AntII_l);
	AntIII_IV=mean(of AntIII_IV_r, AntIII_IV_l);
	AntV=mean(of AntV_r, AntV_l);
	VIb=mean(of VIb_r, VIb_l);
	pt=mean(of pt_r, pt_l);
	AntIIIBW=mean (of AntIIIBW_r, AntIIIBW_l);
	NsetaeAntIII_IV=mean(of NsetaeAntIII_IV_r, NsetaeAntIII_IV_l);
	NsetaeAntV=mean(of NsetaeAntV_r, NsetaeAntV_l);
	Ant=mean(of ant_r, ant_l);
	pt_cauda=pt/cauda;
	HindTibia_pt=HindTibia/pt;
	Ant_BL=Ant/BL;
	urs_htII=urs/htII;
	pt_VIb=pt/VIb;
	pt_siph=pt/siph;
	cauda_urs=cauda/urs;
	siph_BL=siph/BL;
	siph_siphBW=siph/siphBW;
	siph_cauda=siph/cauda;
	urs_VIb=urs/VIb;
	siph_HindFemur=siph/HindFemur;
	cauda_caudaBW=cauda/caudaBW;
	drop htII_r htII_l HindFemur_r HindFemur_l HindTibia_r HindTibia_l HindTibiaW_r
	HindTibiaW_l siph_r siphDW_r siphBW_r siph_l siphDW_l siphBW_l AntI_l AntII_l AntIII_l AntIIIBW_l
	AntIV_l AntV_l VIb_l pt_l Ant_l AntI_r AntII_r AntIII_r AntIIIBW_r AntIII_IV_r AntIII_IV_l
	AntIV_r AntV_r VIb_r pt_r Ant_r setaeAntIII_IV_r NsetaeAntIII_IV_l NsetaeAntV_r NsetaeAntV_l
	;
run;
```
### Data analysis
#### synthesis table (Supp. Table 1): computation of mean, min, max, by species
```
data morphometry_synth;
  set morpho2;
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
