# deBruijn
The project takes target database (.fasta) file as input and generates a de Bruijn decoy (.fasta) file and a target+decoy concatenated (.fasta) file.



de Bruijn decoy parameters
--------------------------
### Options:

#### Required argument:

-  **--input <.fasta>** 

#### Optinal arguments:

-  **[--output <.fasta>]**

-  **[--k <integer value greater than 0>] (default: 2)**

-  **[--enzyme 0/1/2] (0=non specific, 1=trypsin (default), 2=chymotrypsin)**

- **[--compositoin 0/1/2] (0=not adaptive, target composition not considered, 1=not adaptive, target composition considered, 2=adaptive, target composition consodered (default))**


Running de Bruijn decoy generator
---------------------------------

- Download the ZIP file and unzip deBruijn-master.zip

- Open the project in an IDE (e.g., IntelliJ or NetBeans), build the project

- Run the project from the IDE or

### Command line

- Build: ant -f "DIRECTORY OF THE DEBRUIJN-MASTER" -Dnb.internal.action.name=rebuild clean jar

- Run: java -jar "DIRECTORY OF THE DEBRUIJN-MASTER/dist/deBruijn.jar --input <*.fasta> [optional arguments]"

