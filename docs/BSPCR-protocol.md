# Bisulfite PCR protocol
## Background
Bisulfite conversion of DNA is a method of deaminating unmethylated cytosines into uracil, leaving only 5-methylcytosine (5mC) and its other forms intact. This technique is used in the analysis of methylation status of DNA. By performing the bisulfite conversion, sequencing the converted DNA, and comparing the converted sequence to an unconverted reference, we can identify which cytosines are methylated, as they will be read as cytosine while unconverted cytosines will be read as thymine. 

Note that the conversion process causes both strands of DNA to no longer be complementary to each other. This is due to unmethylated cytosines getting converted to uracil (and subsequently thymine after PCR), resulting in G-T mismatches across both strands. Therefore, bisulfite-converted DNA is single-stranded and is significantly less stable. Bisulfite-converted DNA also has the tendency to form complex secondary structures due to its single-strandedness. Finally, bisulfite-conversion is damaging to DNA due to the high alkalinity and temperatures required in the reaction. As a result, bisulfite-converted DNA will often be fragmented.

Due to the unique characteristics of bisulfite-converted DNA, PCRs on bisulfite-treated DNA require additional considerations: 
* Primers for PCRs on bisulfite-treated DNA must be designed to target bisulfite-treated DNA. Ideally, primers should not have any CG dinucleotides in them, as these CGs may be either be a C or a T depending on their methylation status. If CGs are present in the primers, they should be placed on the 5’ end of the primer. Conventional primer design programs cannot be used to design these primers. See the ‘Bisulfite PCR Primers Design using Bisearch’ protocol for more details. 

* Since both strands of DNA are non-complementary and unique after bisulfite-conversion, any primers that are designed will only amplify off one strand. Therefore, to completely assess the methylation status of any region, at least two sets of primers are required. 

* The complex secondary structures formed by bisulfite-converted DNA often mean that PCRs will require the use of PCR additives and buffers that disrupt base-pairing and secondary structure formation. In my experience, both betaine and DMSO are both needed. For Phusion U, using the Phusion GC buffer significantly increases the success rate of the BSPCR, as the GC buffer is also designed to disrupt base-pairing and secondary structure formation. 

* Because of bisulfite-converted DNA is fragmented, there is an upper limit to the size of amplicon that can be produced by PCR. Generally, good yields will be obtained from PCRs of <300 bp fragments, with worse and worse yields up to ~600 bp. 

* Fragmentation as well as incomplete bisulfite conversion of priming sites reduce the number of amplifiable copies of DNA in your sample. Therefore, expect that yields of PCRs on bisulfite-converted DNA to be significantly less than standard PCRs, and will often require additional cycles to produce a usable amount of product. 

* The deamination of cytosine produces uracil. Copies of the region with thymine are produced during the PCR. Therefore, the polymerase used must be able to amplify off uracil-containing templates. Refer to the manufacturer’s product information sheets to identify if your polymerases are suitable for amplification off bisulfite-converted DNA.

* PCRs on bisulfite-converted DNA are often more sensitive to changes to annealing temperature than conventional PCRs. During the optimization of primer pairs for bisulfite-converted DNA, I would recommend optimizing with 1°C increments around the calculated annealing temperature of the primer pair. 

## Equipment and Reagents
### Reagents
| Item | Cat. # | Concentration | Notes |
| ---- | ------ | ------------- | ----- |
| Phusion U Polymerase | F555S/L | | |
| Invitrogen UltraPure Agarose | 16500-500 | N/A | |
| Invitrogen SYBR Safe DNA Gel Stain | S33102 | 20000x | |
| DMSO | N/A | 100% | |
| dNTPs | N/A | 2 mM | |
| Betaine | N/A | 5 M | |
| Primers | N/A | 5 uM | Forward and reverse primers |
| Template DNA | N/A | >0.25 ng/uL | |

### Equipment/Consumables
| Equipment | Consumables |
| --------- | ----------- |
| Thermocyclers | PCR tubes/strips |
| Micropipettes | Pipette tips |
| Microcentrifuge | Microcentrifuge tubes |
| Table-top vortex | |

## Protocol 
1. Prepare a master-mix for your PCR. Below is a recommended starting recipe for a 10 uL reaction volume BSPCR. 

    #### Table 1: BSPCR master-mix
    | Reagent                      | Volume (uL)  | Final Concentration (10 uL Rxn) |
    | ---------------------------- | ------------ | ------------------------------- |
    | Phusion GC Buffer, 5X        | 2            | 1X                              |
    | Betaine, 5M                  | 2            | 1M                              |
    | dNTPs, 2mM                   | 1            | 0.8 mM                          |
    | Nuclease-free water          | 0.6          | \-                              |
    | DMSO                         | 0.3          | 3%                              |
    | FW & REV Primer, 5 µM        | 1 (for each) | 0.5 µM                          |
    | Phusion U Polymerase, 2 U/uL | 0.1          | 0.2 U                           |
    
    When preparing the master mix, combine all components other than the polymerase and vortex well. Afterwards, add the polymerase and mix well via flicking and inversion. Vortexing is also possible but keep the pulses as small as possible (e.g. 0.5 seconds). Centrifuge after mixing. 
2. Deliver the appropriate volume of master-mix for each reaction you are preparing. If the master-mix contains the primers, then 8 uL of master-mix should be aliquoted. If not, 6 uL of master-mix should be aliquoted, assuming that you are adding 1 uL of 5 uM forward and reverse primers, respectively.
3. Add primers then the template DNA (2 uL, dilute if necessary) into each aliquot. 
4. Mix the PCR reactions well via vortexing or flicking and inverting the tubes, then centrifuge.
    * Vortexing briefly seems to produce more consistent results than flicking and inverting the tubes. 
5. Load the tubes into a thermocycler, and use the following program: 
    #### Table 3: Thermocycling protocol for Phusion U 
    | Temperature | Duration  |
    | ----------- | --------- |
    | 98°C        | 30s       |
    | 98°C        | 10s (x40) |
    | X°C         | 30s (x40) |
    | 72°C        | 30s (x40) |
    | 72°C        | 10 min    |
6. To confirm the success of the PCR, perform agarose gel electrophoresis. Prepare a 2% gel. Load 3 uL of PCR product mixed with 3 uL of loading dye into each well. Load 1 uL of ladder mixed with 3 uL of loading dye and 3 uL of water. Run the gel at 120 V for 40 minutes. 
7. For PCR optimization, prepare 6-9 replicates of each PCR and test them in 1 C increments in the Quantstudio 5. These PCR sets are extremely sensitivie to the annealing temperature. 

## References
1. Thermoscientific. (n.d.). Phusion U Hot Start. [Link.](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012916_PhusionUHotStart_DNAPolymerase_F555L_UG.pdf)
	
2. Darst, R. P., Pardo, C. E., Ai, L., Brown, K. D., & Kladde, M. P. (2010). Bisulfite sequencing of DNA. Current Protocols in Molecular Biology / Edited by Frederick M. Ausubel... [et Al.], Chapter 7, Unit 7.9.1–17.