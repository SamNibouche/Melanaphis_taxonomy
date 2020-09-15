## Data analysis
### Data
Morphometric data are available at [DOI:10.18167/DVN1/PDPDS4](http://dx.doi.org/10.18167/DVN1/PDPDS4)
### Import data
```
data WORK.MORPHOMETRY    ;
	infile 'D:\Mes Données\Etudes\SCYLV\diversité vecteur\morphométrie & taxo\routines stats et data pour dépôt\morphometry_alates_raw.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
	informat specimen $15. ;
	informat species $8. ;
	format specimen $15. ;
	format species $8. ;
	input
		specimen $
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
		NsetaeAntIII_IV_l
		NsetaeAntV_r
		NsetaeAntV_l
		NRhinAntIII_r
		NRhinAntIV_r
		NRhinAntV_r
		NRhinAntIII_l
		NRhinAntIV_l
		NRhinAntV_l
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
	NRhinAntIII=mean(of NRhinAntIII_r, NRhinAntIII_l);
	NRhinAntIV=mean(of NRhinAntIV_r, NRhinAntIV_l);
	NRhinAntV=mean(of NRhinAntV_r, NRhinAntV_l);
	Ant=mean(of ant_r, ant_l);
	Ant_BL=Ant/BL;
	urs_htII=urs/htII;
	pt_VIb=pt/VIb;
	pt_cauda=pt/cauda;
	pt_siph=pt/siph;
	HindTibia_pt=HindTibia/pt;
	cauda_urs=cauda/urs;
	urs_VIb=urs/VIb;
	siph_BL=siph/BL;
	siph_HindFemur=siph/HindFemur;
	siph_cauda=siph/cauda;
	cauda_caudaBW=cauda/caudaBW;
	drop htII_r htII_l HindFemur_r HindFemur_l HindTibia_r HindTibia_l HindTibiaW_r
	HindTibiaW_l siph_r siphDW_r siphBW_r siph_l siphDW_l siphBW_l AntI_l AntII_l AntIII_l AntIIIBW_l	
	AntIV_l	AntV_l VIb_l pt_l Ant_l AntI_r AntII_r AntIII_r AntIIIBW_r AntIII_IV_r AntIII_IV_l	
	AntIV_r	AntV_r VIb_r pt_r Ant_r NsetaeAntIII_IV_r NsetaeAntIII_IV_l NsetaeAntV_r	
	NsetaeAntV_l NRhinAntIII_r NRhinAntIV_r NRhinAntV_r NRhinAntIII_l NRhinAntIV_l NRhinAntV_l NRhinAntIV_r
	;
run;

```
### Data analysis
#### synthesis table (Supp. Table 2): computation of mean, min, max, by species
```
data morphometry_synth;
  set morpho2;
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

