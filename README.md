# GeneScoring

## Description 

GeneScoring is a script that takes a .gmt formatted file and a STRING.txt database and calculates a score of each gene
found in the .gmt formatted file. The score is calculated by taking the average number of connections that the gene
has in n_trials subnetworks. A .sif file is then written to the specified directory wherein the genes containing the
highest scores in each loci make up the subnetwork. Visualization is done in Cytoscape.

## Workflow

### Arguments 

Here is a list of required and optional arguments that can be accessed by 
typing **python "GeneScoring.py" -h**: 

```text
usage: GeneScoring.py [-h] [--nt NT] [--gmt GMT] [--sdb SDB] [--sfd SFD]

Calculate average gene scores and visualize networks of the genes with the highest scores from each loci

optional arguments:
  -h, --help            show this help message and exit
  --nt NT, -n_trials NT
                        number of trials to run
  --gmt GMT, -gmt_formated_file GMT
                        experimental loci file path (default "Input.gmt.txt")
  --sdb SDB, -string_database SDB
                        gene interaction file path (default "STRING.txt")
  --sfd SFD, -SIF_file_directory SFD
                        file path of new sif file (default "subnetwork.sif")
```
### Inputs 

Requires one file STRING.txt database, a .gmt format text file, an integer
for the number of trials you want to run (default 5000), and the file path for new SIF file. 

Input .gmt format file example: 

```text
Fanconi anemia locus 0	Locus for PALB2	NUPR1	CTB-134H23.2	SLC5A11	KIAA0556	CD19	SH2B1	CCDC101	GTF3C1	IL27	ARHGAP17	ERN2	DCTN5	NSMCE1	AQP8	RABEP2	XPO6	ATP2A1	CHP2	BOLA2	KDM8	EIF3C	ATXN2L	LAT	ZKSCAN2	SULT1A1	HS3ST4	EIF3CL	TUFM	NPIPL1	SNX29P2	IL21R	PRKCB	SPNS1	TNRC6A	CACNG3	PLK1	RBBP6	NFATC2IP	APOBR	IL4R	PALB2	SULT1A2	CTD-3203P2.2	GSG1L	SBK1	LCMT1
Fanconi anemia locus 1	Locus for FANCF	CSTF3	FBXO3	SLC17A6	CCDC73	CAPRIN1	RCN1	BDNF	METTL15	CCDC34	EIF3M	LUZP2	BBOX1	CAT	PRRG4	SLC5A12	QSER1	AC103801.2	TCP11L1	SVIP	CD59	NAT10	C11orf91	KCNA4	FIBIN	ARL14EP	ABTB2	LMO2	ELP4	MUC15	DNAJC24	GAS2	LGR4	RP11-17A1.2	MPPED2	KIF18A	FANCF	FSHB	HIPK3	PAX6	DEPDC7	IMMP1L	KIAA1549L	WT1	DCDC5	AC132216.1	ANO5	ANO3	ELF5	EHF	LIN7C
Fanconi anemia locus 2	Locus for RAD51C, BRIP1	CD79B	CACNG1	TANC2	SMG8	RP11-15E18.4	TEX2	YPEL2	RGS9	C17orf72	STRADA	DDX42	TACO1	ICAM2	APOH	PRKCA	FTSJ3	TBX4	DCAF7	GH1	GDPD1	CTD-2510F5.6	METTL2A	MRC2	MAP3K3	PRR11	MED13	C17orf64	TBX2	POLG2	SMURF2	AXIN2	CEP95	INTS2	RAD51C	PPM1E	CA4	CEP112	SMARCD2	C17orf82	USP32	KCNH6	CACNG4	CSH1	RPS6KB1	CTD-2535L24.2	RNFT1	BCAS3	LIMD2	NACA2	RP11-51F16.8	DDX5	APPBP2	SKA2	TRIM37	SCN4A	PTRH2	DHX40	RP11-178C3.1	CCDC47	GNA13	GH2	CSH2	CYB561	HEATR6	VMP1	PSMC5	CSHL1	EFCAB3	TUBD1	CACNG5	BRIP1	PPM1D	AC005544.1	LRRC37A3	CLTC	ERN1	MARCH10	TLK2
Fanconi anemia locus 3	Locus for FANCC	TMOD1	MSANTD3-TMEFF1	TEX10	HIATL2	LPPR1	MRPL50	CDC14B	FOXE1	ZNF510	BAAT	COL15A1	CORO2A	HABP4	INVS	SLC35D2	C9orf156	HSD17B3	ZNF189	GABBR2	C9orf174	FAM22G	ZNF782	ANP32B	XPA	DKFZP434H0512	TSTD2	SEC61B	NR4A3	ANKS6	PTCH1	TMEFF1	MURC	CTSL2	AAED1	STX17	GALNT12	ERP44	TGFBR1	TDRD7	TBC1D2	FANCC	HEMGN	ALG2	MSANTD3	TRIM14	NANS	NCBP1	KRT8P11	ERCC6L2	ZNF367
Fanconi anemia locus 4	Locus for FANCA	DBNDD1	AC133919.6	RP11-356C4.2	MC1R	SPIRE2	C16orf3	CENPBD1	FANCA	DEF8	GAS8	PRDM7	TCF25
Fanconi anemia locus 5	Locus for UBE2T	PPFIA4	SNRPE	LGR6	ZC3H11A	TMEM183A	PPP1R12B	SOX13	GOLT1A	DSTYK	NUAK2	NFASC	LRRN2	ADIPOR1	PPP1R15B	UBE2T	CHIT1	KLHL12	MDM4	ZBED6	RBBP5	RABIF	BTG2	MYOG	KISS1	REN	CDK18	PIK3C2B	KDM5B	CNTN2	CYB5R1	LAX1	ADORA1	ELK4	TMCC2	AC119673.1	KLHDC8A	RP11-480I12.4	MFSD4	MYBPH	ETNK2	OPTC	SYT2	FMOD	PLEKHA6	CHI3L1	PRELP	LEMD1	TMEM81	ATP2B4	SLC45A3
Fanconi anemia locus 6	Locus for FANCD2	VGLL4	TAMM41	ZFYVE20	EAF1	IQSEC1	SLC6A6	VHL	FANCD2	TMEM43	C3orf20	TATDN2	PPARG	LSM3	CCDC174	NUP210	SH3BP5	WNT7A	TMEM40	HRH1	IRAK2	SEC13	TSEN2	ATP2B2	SLC6A11	CHCHD4	COLQ	MKRN2-AS1	BRK1	CAPN7	AC034193.1	FGD5	TPRXL	XPC	GHRL	SYN2	ATG7	MKRN2	SLC6A1	NR2C2	BTD	FANCD2OS	MRPS25	HDAC11	HACL1	RPL32	RAF1	METTL6	TIMP4	CAND2	FBLN2
Fanconi anemia locus 7	Locus for FANCE	KCNK5	MTCH1	PI16	BRPF3	TMEM217	RPL10A	GLP1R	ARMC12	STK38	KIF6	CLPSL2	CCDC167	DAAM2	PPIL1	MDGA1	MAPK13	KCNK17	RNF8	MAPK14	ETV7	SRSF3	PXT1	COX6A1P2	TEAD3	C6orf89	FTSJD2	SAYSD1	SRPK1	KCNK16	TULP1	ZFAND3	CPNE5	PPARD	GLO1	FKBP5	PNPLA1	TBC1D22B	SLC26A8	KCTD20	C6orf222	BTBD9	PIM1	FANCE	LHFPL5	CDKN1A	CLPSL1	DNAH8	RAB44	FGD2	CLPS
Fanconi anemia locus 8	Locus for BRCA2	PROSER1	MRPS31	UFM1	AKAP11	N4BP2L2	TNFSF11	FOXO1	CCDC169	KBTBD7	PDS5B	EXOSC8	FREM2	MAB21L1	VWA8	TRPC4	NHLRC3	POSTN	RFC3	N4BP2L1	STARD13	FAM216B	BRCA2	DCLK1	MTRF1	LHFP	ELF1	KL	CSNK1A1L	SOHLH2	DGKH	NAA16	NBEA	SMAD9	FAM48A	STOML3	SPG20	CCDC169-SOHLH2	SLC25A15	RFXAP	SERTM1	COG6	RGCC	RP11-298P3.4	SPG20OS	KBTBD6	SUGT1P3	ALG5	AL133318.1	WBP4	CCNA1
Fanconi anemia locus 9	Locus for ERCC4	GPR139	C16orf88	RPS15A	BFAR	CTD-2349B8.1	NDE1	PLA2G10	RP11-1212A22.2	SYT17	GPRC5B	FOPNL	ITPRIPL2	PDXDC1	UMOD	GDE1	RP11-467M13.1	COQ7	RP11-719K4.1	TMC7	C16orf62	KIAA0430	SMG1	NOMO2	NTAN1	ARL6IP1	PARN	CLEC19A	RRN3	XYLT1	IQCK	NOMO1	ABCC1	C16orf45	CCP110	NPIP	TMC5	RP11-1035H13.3	MKL2	RP11-719K4.2	GP2	NPIPP1	PDILT	AC092291.2	SHISA9	MYH11	NOMO3	ABCC6	ERCC4
Fanconi anemia locus 10	Locus for FANCI	POLG	AC016251.1	FANCI	GABARAPL3	FES	C15orf38	VPS33B	NR2F2	SLCO3A1	IQGAP1	FAM174B	C15orf32	CRTC3	HDDC3	MESP2	GDPGP1	ANPEP	MCTP2	MAN2A2	TTLL13	SPATA8	BLM	ST8SIA2	TICRR	C15orf38-AP3S2	UNC45A	RCCD1	NGRN	SEMA4B	CIB1	AC112693.2	SV2B	ARRDC4	RGMA	PRC1	MESP1	KIF7	CHD2	IDH2	PLIN1	RP11-697E2.6	PEX11A	WDR93	RHCG	ZNF710	ZNF774	RP11-82I10.1	FURIN	AC087477.1	AP3S2
Fanconi anemia locus 11	Locus for SLX4	FAM86A	USP7	ROGDI	ALG1	CDIP1	UBN1	RP11-297M9.1	SEPT12	CORO7	NUDT16L1	C16orf71	PAM16	TMEM186	TRAP1	METTL22	SRL	RBFOX1	VASN	GLYR1	CARHSP1	GRIN2A	C16orf89	ATF7IP2	NAGPA	DNASE1	FAM100A	TEKT5	SEC14L5	AC012676.1	PMM2	GLIS2	C16orf96	ZNF500	C16orf72	CORO7-PAM16	DNAJA3	CREBBP	NMRAL1	ADCY9	NLRC3	HMOX2	MGRN1	EMP2	RP11-127I20.4	ABAT	SLX4	PPL	TMEM114	TFAP4	ANKS3
```
Input STRING.txt database example: 
```
ARF5	DVL2	0.166000
ARF5	DYRK4	0.166000
ARF5	PPP5C	0.254968
ARF5	MAP4K5	0.157276
ARF5	RALBP1	0.156000
ARF5	PKP2	0.160210
ARF5	ACAP1	0.328000
ARF5	MAP2K5	0.242000
ARF5	MYO15A	0.272395
ARF5	MAPK13	0.190000
ARF5	STX1B	0.263160
ARF5	MAPK12	0.190000
ARF5	MAPK1	0.190000
ARF5	MYH9	0.252822
ARF5	SOS2	0.199000
ARF5	EIF5	0.214358
ARF5	PABPC1L	0.196000
ARF5	CORO1A	0.163000
ARF5	PARD6A	0.194000
ARF5	PDIA2	0.157397
...
```
### Command

To run this program in the command line interface type: 
```text
python "GeneScoring.py" --nt 5000 --gmt "Input.gmt.txt" --sdb "STRING.txt" --sfd "subnetwork.sif"
```
You can replace --gmt and --sdb with your own desired input files as well as change 
the number of trials (--nt) and the file name (--sfd). The file will be automatically put in 
the current working directory unless otherwise specified.

To run this program interactively type: 


```python
VisualizeGeneScoring(5000, "Input.gmt.txt", "STRING.txt", "subnetwork.sif")
```
If the input files are located in a different directory then you can put their respective file paths. 

### Output 

Example of SIF file formatting:

```text
PLK1	pp	RAD51C	CDK18	DCLK1	ERCC4	TRAP1
LGR4	pp	NCBP1	MAPK14	TRAP1
RAD51C	pp	PLK1	CDK18	DCLK1	ERCC4
NCBP1	pp	LGR4	MAPK14
FANCA	pp	ERCC4
CDK18	pp	PLK1	RAD51C	IRAK2	MAPK14	DCLK1	ERCC4	TRAP1
IRAK2	pp	CDK18	MAPK14	CIB1
MAPK14	pp	LGR4	NCBP1	CDK18	IRAK2	DCLK1	ERCC4	CIB1	TRAP1
DCLK1	pp	PLK1	RAD51C	CDK18	MAPK14	ERCC4	CIB1
ERCC4	pp	PLK1	RAD51C	FANCA	CDK18	MAPK14	DCLK1
CIB1	pp	IRAK2	MAPK14	DCLK1
TRAP1	pp	PLK1	LGR4	CDK18	MAPK14
```
### Visualization 
Visualization was done in Cytoscape. The following is the visualization of the example SIF file: 

![Image of subnetwork example](https://github.com/acolorado1/Module3Day3/blob/f1d3c6f64d2ff337f12c23b13ee7be9f235d0b07/Example1Subnetwork.png)

## Installation and Dependencies
You must have Python 3 installed. Any Python 3 version should work but it was written in Python 3.9 using a Windows-based 
operating system. Packages random, argparse 1.4.0, and the function mean from the statistics package will need to be 
installed. 

## Contact 
Angela Sofia Burkhart Colorado - angelasofia.burkhartcolorado@cuanschutz.edu

Project Link: https://github.com/acolorado1/Module3Day3.git