/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deBruijn;

import java.util.HashMap;
import java.util.Random;
import java.util.ArrayList;
import java.util.Comparator;
//import java.util.concurrent.TimeUnit;

import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
//import org.apache.commons.cli.*;

/**
 *
 * @author jmmoosa
 */
public class deBruijn {
    
    static final double WATER_MASS=18.01;//0565;
    static final double PROTON_MASS=1.0073;

    /**
     * list of Amino Acids (alphabets) in different enzyme conditions
     */
    public static final String NON_SPECIFIC_AA="ACDEFGHIKLMNPQRSTVWY";//20 amino acids
    public static final String ALL_A="ACDEFGHIKLMNPQRSTVWYBJOUXZ";//26 //just a list of all alphabets
    public static final String TRYPSIN_AA="ACDEFGHILMNPQSTVWY";//18, without K and R, 
    public static final String CHYMOTRYPSIN_AA="ACDEGHIKLMNPQRSTV";//17, without W, Y, F
    public static final String [] AA_S={NON_SPECIFIC_AA, TRYPSIN_AA, CHYMOTRYPSIN_AA};
    
    public static final int NO_OF_AA=20;
    public static int NO_OF_ENZYME_AA =18;/* Number of AAs, without the cleavage AAs (18 for default enzyme Trypsin) */

    /**
     *Number of alphabets
     */
    public static final int ALPHABET_SIZE=26;
    

    public static final int TRYPSIN=1;
    public static final int CHYMOTRYPSIN=2;
    public static final int NONSPECIFIC=0;

    
    //To format the decoy template output, not necessary
    private final int maximum_entry_per_line=15;     
    //variable to read/(format) input/(output)
    private static int newline=2;
    //tag to insert before the decoys
    public String decoy_tag = "DeBruijn";
    
    
    // Input arguments
    public static int k=2;/* default */
    public static int ENZYME=TRYPSIN; //default
    
    //flags
    //no_KP_RP=true means there will be no new KP/RP introduced in the decoy, no_KP_RP=false means introduction of KP/RP is possible
    public boolean no_KP_RP=true;

    /*consider_target_composition=true means that deBruijn decoy will be generated according to the target Amino Acid (aa) target_composition, 
    consider_target_composition=false means deBruijn decoy will be generated using random aa*/
    public static boolean consider_target_composition=true;

    public static boolean isAdaptive =false;
    
    
    String inputfilepathName, inputfilename, outputfilename;
    
    //class variables
    private readFASTA target_fasta_reader;//, decoy_fasta_reader;
    private writeFASTA decoy_fasta_writer, concat_fasta_writer;
    
    //private read_writeFASTA target_fasta_reader_writer, decoy_fasta_reader_writer;
    private write_DeBruijn de_bruijn_writer;//to write the template
    
    private int[] k_index;
    private int[] k_sub;
   
    
    public HashMap<Character, Double> target_composition, target_composition_wo_fixedAA, decoy_composition, adaptive_composition, normalized_composition;
    public static ArrayList<AminoAcid> cummulative_sum_prob_AA; //stores the cummulative probability of the amino acids, sorted in increasing order
    
    public static final char DUMMY_KP_RP='_'; //dummy chcracter to represent KP/RP
    public static final char DUMMY='B'; //dummy peptide
    
    
    
    HashMap<String, K_MER> de_bruijn_decoy_random, de_bruijn_decoy_input, de_bruijn;
    public static HashMap<String, AminoAcid> AminoAcids; //Key is the single letter Amino Acid code
    
    private static final Logger LOGGER = Logger.getLogger(deBruijn.class.getName());

    /**
     * Name: private void init_list_AA()
     * Purpose: To initialize the list of Amino Acid names and masses (the masses are only needed to generate statistics)
     * Input: None
     * Return: return HashMap<String, AminoAcid>
     * Effects: the public HashMap variable AminoAcids assigned the returned value
    */
    private HashMap<String, AminoAcid> init_list_AA()
    {
        
        AminoAcid aa;     
        
        HashMap<String, AminoAcid> AAs=new HashMap<>();                
        //might need to update the masses
        aa=new AminoAcid("", "", 0.0, 0.0);//dummy
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("G", "Gly", 57.02146, 57.052, 30.03438);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("A", "Ala", 71.03711, 71.079, 44.05003);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("S", "Ser", 87.03203, 87.078, 60.04494);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("P", "Pro", 97.05276, 97.117, 70.06568);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("V", "Val", 99.06841, 99.133, 72.08133);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("T", "Thr", 101.04768, 101.105, 74.06059);//
        AAs.put(aa.single_letter_code, aa);
        
        /*aa=new AminoAcid("C", "Cys", 103.00919, 103.144, 76.0221);////more c available
        AminoAcids.put(aa.single_letter_code, aa);*/
        
        aa=new AminoAcid("C", "Cys", 160.03065, 0.0, 133.0436);////more c available
        AAs.put(aa.single_letter_code, aa);
                
        aa=new AminoAcid("I", "Ile", 113.08406, 113.160, 86.09698);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("L", "Leu", 113.08406, 113.160, 86.09698);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("N", "Asn", 114.04293, 114.104, 87.05584);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("D", "Asp", 115.02694, 115.089, 88.03986);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("Q", "Gln", 128.05858, 128.131, 101.0715);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("K", "Lys", 128.09496, 128.174, 101.1079);////one other immo mass available 84.08136
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("E", "Glu", 129.04259, 129.116, 102.0555 );//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("M", "Met", 131.04049, 131.198, 104.0534);////one other M there
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("H", "His", 137.05891, 137.142, 110.0718);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("F", "Phe", 147.06841, 147.177, 120.0813);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("R", "Arg", 156.10111, 156.188, 129.114);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("Y", "Tyr", 163.06333, 163.170, 136.0762);//
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("W", "Trp", 186.07931, 186.213, 159.0922);//
        AAs.put(aa.single_letter_code, aa);
        
        //x & u, not sure
        
        aa=new AminoAcid("U", "Sec", 150.95363, 150.0379);//U
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("X", "I/L", 113.08406,	113.1594);//X
        AAs.put(aa.single_letter_code, aa);
        //need B, Z
        
        aa=new AminoAcid("B", "Asx", 114.04293, 114.104, 87.05584);//N/D, doing it considering N for now
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("Z", "Gsx", 128.05858, 128.131, 101.0715);//Q/E, doing it considering Q for now
        AAs.put(aa.single_letter_code, aa);
        
        aa=new AminoAcid("O", "Pyl", 255.31, 255.31, 255.31);//O just from wiki
        AAs.put(aa.single_letter_code, aa);
        
        return AAs;
        
    }


    deBruijn (String inputDirectory, String inputfileName, String outputfilename, int k_a)
    {
        //the constructor assignments
        k=k_a;
        this.inputfilepathName=inputDirectory+inputfileName;
        this.outputfilename=outputfilename;


        //load the AminoAcids in this HashMap variable
        AminoAcids = init_list_AA();
        target_composition=new HashMap<>();
        target_composition = init_composition_hashmap(); //to initialize the variable, the hashmap named "target_composition"


        cummulative_sum_prob_AA =new ArrayList<>();

        de_bruijn_decoy_random=new HashMap<>();
        de_bruijn_decoy_input=new HashMap<>();
        de_bruijn=new HashMap<>();

        de_bruijn_writer=new write_DeBruijn(inputDirectory+"e"+ENZYME+"k"+k+"deBruijn_template.txt");

        k_index=new int[k+1];
        k_sub=new int[k+1];
        
        NO_OF_ENZYME_AA = AA_S[ENZYME].length();
        
        //file reader writer initialization
        target_fasta_reader = new readFASTA(this.inputfilepathName);
        //decoy_fasta_reader = new readFASTA(outputfilename);
        decoy_fasta_writer = new writeFASTA(inputDirectory+outputfilename);
        concat_fasta_writer = new writeFASTA(inputDirectory+"concat-MQ-"+outputfilename);/* not sure if needed? maybe add a parameter option? */
        
        //

        adaptive_composition=new HashMap<>();
        decoy_composition=new HashMap<>();        
        target_composition_wo_fixedAA = new HashMap<>();
        normalized_composition = new HashMap<>();

        getLoggerFile();
    }
    
    // to initialize the hashmap named "target_composition", with the AminoAcids, intial target_composition is set to 0.0
    /**
     * Name: private HashMap<Character, Double> init_composition_hashmap()
     * Input: None
     * Return: HashMap<Character, Double>
     * Effect: the returned value is assigned to the variable target_composition
    */
    private HashMap<Character, Double> init_composition_hashmap()
    {
        HashMap<Character, Double> composition=new HashMap<>();

        for(int i=0; i<ALPHABET_SIZE; i++)//initialized for all possible alphabets
        {
            char aa=ALL_A.charAt(i);
            composition.put(aa, 0.0); //initialize compositions to 0.0
        }

        composition.put(DUMMY_KP_RP, 0.0); //to store number of KP/RP

        return composition;
    }

    /**
     * Purpose: To copy a composition from one HashMap to other.
     * @param source: the source composition, from where the composition will be copied
     * @param dest: the destination, to where the composition will be copied
     * @return dest HashMap after the composition is copied
     */
    private HashMap<Character, Double> copy_composition_to_hashmap(HashMap<Character, Double> source, HashMap<Character, Double> dest)
    {     
        for(Character aa : source.keySet())
        {
            double count = source.get(aa);
         
            dest.put(aa, count);
        }
        
        return dest;
    }

    /**
     *
     * @param fasta_reader the reader to be opened, not necessary, as only one line. no need for method?
     */
    private void open_reader(readFASTA fasta_reader)
    {
        fasta_reader.openReader();        
    }

    /**
     *
     * @param fasta_reader the reader to be opened
     * @param fasta_writer the writer to be opened
     * Effect: also opens the concat_fasta_writer to write the concatenated fasta file, and de_bruijn_writer to write the template text file
     */
    private void open_reader_writer(readFASTA fasta_reader, writeFASTA fasta_writer)
    {
        fasta_reader.openReader();
        fasta_writer.openWriter();
        concat_fasta_writer.openWriter();

        de_bruijn_writer.openWriter();
    }

    /**
     *
     * @param fasta_reader reader file to be closed
     * @param fasta_writer writer file to be closed
     */
    private void closeFile(readFASTA fasta_reader, writeFASTA fasta_writer)
    {        
        closeReader(fasta_reader);
        closeWriter(fasta_writer);
    }

    /**
     *
     * @param fasta_reader reader file to be closed
     */
    private void closeReader(readFASTA fasta_reader)
    {
        fasta_reader.closeReader();
    }

    /**
     *
     * @param fasta_writer writer file to be closed
     * Effect: also the concat_fasta_writer and de_bruijn_writer are closed
     */
    private void closeWriter(writeFASTA fasta_writer)
    {
        fasta_writer.closeWriter();
        concat_fasta_writer.closeWriter();//this
        de_bruijn_writer.closeWriter();

    }

    /**
     * Purpose: generates decoy template for all possible kmers, completely randomly not considering input target_composition
     * @param de_bruijn_decoy_random need to return these? lagbe? necessary?
     */
    void generate_decoy_template_random(HashMap<String, K_MER> de_bruijn_decoy_random)//generates decoy template for all possible kmers, not according to input target_composition, but randomly
    {
        char k_mer[]=new char[k+1];
        
        int no_k_plus1mers =(int)Math.pow(NO_OF_AA, k+1);
        
        for(int j = 0; j< no_k_plus1mers; j++)//generate all possible k-mers
        {         
            k_sub=new int[k+1];
            for(int i=0; i<k+1; i++)
            {
                int index=map_2d_index_1d(i, j);
                //System.out.println(i+", "+j+": "+index);
                k_mer[i]=NON_SPECIFIC_AA.charAt(index);
            }
            
            String kmer = new String (k_mer);
            
            char randomAA;

            //generate the replacement AA
            if(!no_KP_RP)
                randomAA=generate_random();//complete random
            else
            {
                randomAA=generate_random(kmer.charAt(k-1), kmer.charAt(k));// random AA, however no KP/RP generated, K/R reserved
            }
            
            de_bruijn_decoy_random.put(kmer, new K_MER(randomAA)); //all possible k-mers included with a replacement, outgoing not considered
            
        }
    }
            
    //according to input target_composition

    /**
     * Purpose: generates decoy template for all possible kmers, randomly considering input target_composition.
     * Goal is to keep the target and decoy composition similar. So, occurrence of each AA in target and decoy should be similar.
     * @param de_bruijn_decoy_input all the kmers present in the input target fasta, generated during composition calculation
     * @param composition composition of the AAs in the target fasta, contains the total number of occurrence of each AA
     * @param cumsum_prob cumulative probability?
     * @return updated composition
     * Effect: outgoing?
     */
    HashMap<Character, Double> generate_decoy_graph_input_random(HashMap<String, K_MER> de_bruijn_decoy_input, HashMap<Character, Double> composition, ArrayList<AminoAcid> cumsum_prob)
    {
        HashMap<Character, Double> normalized = new HashMap<>();
        
        for(String kmer: de_bruijn_decoy_input.keySet())
        {
            char randomAA;
            int outgoing=de_bruijn_decoy_input.get(kmer).outgoing;//number of positions of the target is going to be replaced with the selected AA
            //first two conditions (if and else {if} works same as generate_decoy_template_random() where replacement is a random AA
            //however, uses the list of kmers from input
            
            if(!no_KP_RP)//no_KP_RP=false, means might introduce new P after K/R in the decoy 
                randomAA=generate_random(); //composition is not considered here
            else
            {
                if(!consider_target_composition)
                    randomAA=generate_random(kmer.charAt(k-1), kmer.charAt(k));// so now no KP/RP generated, K/R reserved, however, target_composition is not considered here
                else
                    if(isAdaptive)//adaptive, the probability of generating an AA updated with assignment
                    {

                        randomAA= generate_random_using_adaptive_composition(kmer, kmer.charAt(k-1), kmer.charAt(k), cumsum_prob);

                        //update cumsum adaptive composition here according to how many kmers of the target sequences are replaced
                        // then return it, to use it next
                        double count = composition.get(randomAA);/* number of occurrences left in the decoy for the selected AA */

                        //outgoing beshi hole kom ekta select korle? bepar na, next time select hobe na, emon? naki outgoing ta weight er moton set kora jay?

                        if(count-outgoing>0)
                            composition.put(randomAA, count-outgoing);/* updated number of occurrences left in the decoy for the selected AA */
                        else
                            composition.put(randomAA, 0.0);//approximate calculation. the difference should not be too much, so no need to select new AA
                
                        normalized=calculate_normalized_composition_wo_cleavage_AA(composition);//updated normalized
                        cumsum_prob=calculate_cum_probability(normalized);/* updated cumulative sum */
                    }
                    else {
                        randomAA= generate_random_using_composition(kmer.charAt(k-1), kmer.charAt(k), cumsum_prob);//compositions are generated considering without updating the cumulative probability
                    }
            }
            
            de_bruijn_decoy_input.put(kmer, new K_MER(randomAA, outgoing)); //all kmers from target input, with replacements included
            //outgoing not needed anymore, still better to save it?
        }
        
        return composition;
        
    }



    //////////////////////////////
    //functions to generate decoy AA for the template
    //////////////////////////////

    /**
     * Purpose: randomly select a replacement aa to replace the current one, target_composition is not considered,
     * Description: no new cleavage site (K/R or W/Y/F) inserted, controlled by the String array AA_S,
     * however new P after [KR](trypsin) or [WYF] (chymotrypsin) which can result in non-cleavage sites
     * @return the replacement AA
     */
    char generate_random()
    {
        Random rand = new Random();
        char randomAA;
        int  tAAindex = rand.nextInt(NO_OF_ENZYME_AA);
        randomAA=AA_S[ENZYME].charAt(tAAindex);
        
        return randomAA;
            
    }

    /**
     * Purpose: to check if the replacement site is a cleavage site or not
     * @param k last (k-th) AA of the kmer
     * @return a flag saying if this a cleavage site or not
     */
    boolean cleavage_site(char k)
    {
       return ((k=='R' || k=='K')&&ENZYME==TRYPSIN)||((k=='Y' || k=='W' || k=='F')&&ENZYME==CHYMOTRYPSIN);
    }

    /**
     * Purpose: randomly select a replacement aa to replace the current one, target_composition is not considered
     * Description: no new cleavage site (K/R or W/Y/F) inserted, controlled by the String array AA_S, does not generate any new P after [KR](trypsin) or [WYF] (chymotrypsin)
     * however, does not keep the existing P after [KR](trypsin) or [WYF] (chymotrypsin) fixed
     * @param k_1 second last ((k-1)-th) AA of the kmer
     * @return the replacement AA
     */
    char generate_random(char k_1)
    {   
        Random rand = new Random();
        char randomAA;
        int  tAAindex = rand.nextInt(NO_OF_ENZYME_AA);
        
        randomAA=AA_S[ENZYME].charAt(tAAindex);//if non speciifc next if will never be executed, this randomAA will be returned
        
        if(cleavage_site(k_1))
        {
            while(randomAA=='P')//if the prev aa is [K/R](trypsin) or [W/Y/F] (chymotrypsin), randomly select a replacement aa until it is not 'P', so no new KP/RP introduced
            {
                tAAindex = rand.nextInt(NO_OF_ENZYME_AA);//does not generate a new cleavage site (K/R or W/Y/F)
                randomAA=AA_S[ENZYME].charAt(tAAindex);               

            }
        }
        return randomAA;            
    }

    /**
     * Purpose: randomly select a replacement aa to replace the current one, target_composition is not considered
     * Description: no new cleavage site (K/R or W/Y/F), and X inserted, controlled by the String array AA_S,
     * does not generate any new P after [KR](trypsin) or [WYF] (chymotrypsin)
     * also keeps the current P after [KR](trypsin) or [WYF] (chymotrypsin) fixed
     * @param k_1 second last ((k-1)-th) AA of the kmer
     * @param k last (k-th) AA of the kmer
     * @return the replacement AA
     */
    char generate_random(char k_1, char k)//randomly without considering composition, cleavage site (K/R or W/Y/F) fixed
    {    
        if(cleavage_site(k) || (cleavage_site(k_1) && k=='P') || k == 'X')/* fix KP/RP/K/R/ or W/Y/F/WP/YP/FP and X */
            return k;           
        return generate_random(k_1);//okay with nonspecific, this extra wrapper call is needed only for specific enzymes
    }

    /**
     * Purpose: finds the amino acid with maximum cumsum_prob which is less than "probability"
     * @param probability the selected probability
     * @param cumsum_prob the cumulative probability distribution from which the AA will be selected
     * @return the selected AA
     */
    AminoAcid find_aa(double probability, ArrayList<AminoAcid> cumsum_prob)
    {     
        boolean found=false;
        int index=0;
        
        AminoAcid selected=cumsum_prob.get(index);//select the first one by default        
        AminoAcid prev=selected;
        
        for (int i = 1; i<cumsum_prob.size()-1; i++)//from 1, because 0 selected first
        {
            AminoAcid aa = cumsum_prob.get(i);
            
            if(probability<=aa.probability&&probability>prev.probability)//if the current probability is bigger than the selected probability, then found
            {
                found=true;
                selected=aa;
                break;//if maximum is found, then no need to look further
            }
            prev=aa;
            
        }

        if(!found && selected.probability==0.0)//when?
        {
            char randomAA=generate_random();            
            selected = new AminoAcid(randomAA, 0.0);
        }

        //if "probability">max (comp), max(comp) aa will be selected??
        return selected;
    }
    
    private AminoAcid get_random_AA(ArrayList<AminoAcid> cumsum_prob)
    {
        Random rand = new Random();
        //Random rand=new Random();
        double random_prob = rand.nextDouble();        
        AminoAcid selected=find_aa(random_prob, cumsum_prob); // select an Amino Acid that's target_composition is in the range//nonspecific, this is returned
        
        return selected;
        
    }

    /**
     * Purpose: randomly select a replacement aa to replace the current one, using the target_composition
     * Description: no new cleavage site (K/R or W/Y/F), and X inserted, controlled by the String array AA_S,
     * does not generate any new P after [KR](trypsin) or [WYF] (chymotrypsin)
     * also keeps the current P after [KR](trypsin) or [WYF] (chymotrypsin) fixed? rakhe?
     * @param k_1 second last ((k-1)-th) AA of the kmer
     * @param k last (k-th) AA of the kmer
     * @param cumsum_prob
     * @return the replacement AA
     */
    char generate_random_using_composition(char k_1, char k, ArrayList<AminoAcid> cumsum_prob)//does not generate any new P after [KR](trypsin) or [WYF] (chymotrypsin), fixed cleavage (K/R or W/Y/F), according to target composition
    {
        if(cleavage_site(k) || (cleavage_site(k_1) && k=='P') || k == 'X')
            return k;
        
        AminoAcid selected=get_random_AA(cumsum_prob); // select an Amino Acid that's target_composition is in the range//nonspecific, this is returned
        
        if(cleavage_site(k_1))
        {
            while(selected.aa=='P')//not going to introduce any new P after K/R
            {
                selected=get_random_AA(cumsum_prob);      
            }
        }

        return selected.aa;//return the amino acid which corresponds to the desired probability
            
    }

    /** same as the previous method??
     * Purpose: randomly select a replacement aa to replace the current one, using the target_composition? where adaptive applied?
     * Description: no new cleavage site (K/R or W/Y/F), and X inserted, controlled by the String array AA_S,
     * does not generate any new P after [KR](trypsin) or [WYF] (chymotrypsin)
     * also keeps the current P after [KR](trypsin) or [WYF] (chymotrypsin) fixed? rakhe?
     * @param k_1 second last ((k-1)-th) AA of the kmer
     * @param k last (k-th) AA of the kmer
     * @param cumsum_prob
     * @return the replacement AA
     */
    char generate_random_using_adaptive_composition(String kmer, char k_1, char k, ArrayList<AminoAcid> cumsum_prob)//does not generate any new P after [KR](trypsin) or [WYF] (chymotrypsin), fixed cleavage (K/R or W/Y/F), according to target composition
    {
        if(cleavage_site(k) || (cleavage_site(k_1) && k=='P') || k == 'X')
            return k;
        
        AminoAcid selected=get_random_AA(cumsum_prob);
        
        if(cleavage_site(k_1))
        {
            while(selected.aa=='P')//not going to introduce any new P after K/R
            {
                selected=get_random_AA(cumsum_prob);
            }
        }
        return selected.aa;//return the amino acid which corresponds to the desired probability            
    }
    

    //needed for preserving KP/RP, this will be incorporated in the template too

    //not necessary anymore template preserves the pattern,

    /**
     * Purpose: while generating decoy check again if we want to mutate this kmer
     * @param aa k-th kemr
     * @return if the AA should be mutated by the replacement AA
     */
    boolean check_pattern(char aa)//returns true if pattern is not a cleavage site (K/R or W/Y/F)
    {
        boolean isMutationPattern=true;
        if(cleavage_site(aa))
            isMutationPattern=false;
        return isMutationPattern;
    }

    /**
     * Purpose: while generating decoy check again if we want to mutate this kmer
     * @param prev_aa k-1 th kmer
     * @param aa k-th kemr
     * @return if the AA should be mutated by the replacement AA
     */
    boolean check_pattern(char prev_aa, char aa)//returns true if pattern is not a cleavage site (K/R or W/Y/F), if P after [KR](trypsin) or [WYF] (chymotrypsin) return true
    {
        boolean isMutationPattern=true;
        
        if(cleavage_site(aa))
            isMutationPattern=false;
        if(cleavage_site(prev_aa)&& aa=='P')
            isMutationPattern=false;
        
        return isMutationPattern;
    }

    /**
     * Purpose: apply template to the kmer to get the replacement AA
     * @param kmer
     * @return the replacement AA
     */
    char apply_template(String kmer)
    {
        char a=kmer.charAt(k);
        
        if(de_bruijn.containsKey(kmer))
        {
            a=de_bruijn.get(kmer).decoy;
        }
        return a;        
    }

    /**
     * Purpose: wrapper to generate different decoys, not necessary anymore. we only have one decoy
     * @param protein target protein
     * @param newline_positions newline positions in target, to keep target and decoy formatting same
     */
    void generate_decoy(String protein, ArrayList <Integer> newline_positions)
    {
        generate_decoy_de_Bruijn(protein, newline_positions);
    }

    /**
     * Purpose: needed for Shifted Reversal ? lagbe r?
     * @param protein source protein to be modified
     * @param index the index to swap with its previous index
     * @return the modified protein
     */
    String swap (String protein, int index)
    {
        char aa = protein.charAt(index);
        
        if (index>0)
        {
            char[] protein_chars = protein.toCharArray();
            protein_chars[index] = protein_chars[index-1];
            protein_chars[index-1] = aa;
            protein = String.valueOf(protein_chars);
        }
        return protein;
    }


    /**
     * Purpose: to generate deBruijn decoy, and write it
     * @param protein target protein
     * @param newline_positions  newline positions in target, to keep target and decoy formatting same
     */
    void generate_decoy_de_Bruijn(String protein, ArrayList <Integer> newline_positions)
    {
        String kmer="";
        int length=protein.length();
        int j=0;
        int nl=newline_positions.get(j);//needed for formatting only
        
        for(int i=0; i<length-(k+1)+1; i++)
        {
            if(i==nl)
            {
                j++;
                nl=newline_positions.get(j);
                decoy_fasta_writer.write(System.lineSeparator());//this is done only to make the decoy and target look similar in formatting
                concat_fasta_writer.write(System.lineSeparator());//this is done only to make the decoy and target look similar in formatting
            }
            
            kmer=protein.substring(i, i+k+1);
            
            char aa=kmer.charAt(k);
            char prev_aa=kmer.charAt(k-1);
            char next_aa=DUMMY;
            
            if(i<length-(k+1))/* check if last aa */
                next_aa=protein.charAt(i+k+1);
            
            //boolean can_mutate=check_pattern(prev_aa, aa, next_aa);//check korte hobe current er shathe next
            
            boolean can_mutate;
            
            //not necessary anymore template preserves this
            if(!no_KP_RP)
                can_mutate=check_pattern(aa);//if cleavage site (K/R or W/Y/F) then true
            else
                can_mutate=check_pattern(prev_aa, aa);//check if cleavage site (K/R or W/Y/F) or P after [KR](trypsin) or [WYF] (chymotrypsin), if found can_mutate is false

            //char randomAA=kmer.charAt(k);
            if(can_mutate)
            {
                char decoy=apply_template(kmer);
                aa=decoy;
            }
            
            decoy_fasta_writer.write(aa);
            concat_fasta_writer.write(aa);

        }
                
    }

    /**
     * Purpose: among all possible k-mers which kmer is the current one, for each kmers, for each k value which AA should this position get
     * each position has 20 choices
     * @param i one of the 2-d indices, maximum value k
     * @param j one of the 2-d indices, maximum value no_k_plus1mers
     * @return the converted 1-d index, pointing to an AA
     */
    int map_2d_index_1d(int i, int j)
    {
        if(i==0)
        {
            k_index[i]=(int)(j/(int)Math.pow(NO_OF_AA, k-i));
            k_sub[i]=k_index[i]*(int)Math.pow(NO_OF_AA, k-i);
        }
        else
        {
            int sub=k_sub[i-1];
            int n=j-sub;
            k_index[i]=(n/(int)Math.pow(NO_OF_AA, k-i));
            k_sub[i]=sub+(k_index[i]*(int)Math.pow(NO_OF_AA, k-i));

        }
        return k_index[i];
    }

    /**
     * Purpose: read each protein from the target fasta file, generate decoy by applying the pre-generated template, write the decoy fasta
     */
    private void read_generate_write()
    {
        String annotation = "";
        String line, protein="";
        boolean validity=false;
        ArrayList <Integer> newline_positions=new ArrayList<>();

        while((line=target_fasta_reader.readline()) != null)
        {
            if (line.length( ) > 0)//check if empty read
                if(line.charAt(0)!='>')//check if an annotation is found, i.e, a new protein is started
                {
                    validity=true;
                    protein+=line;//concat the whole protein in a single variable
                    newline_positions.add(protein.length()-k);//to keep the same format newline/linebreak positions are saved//-newline (see if needed)
                    //write the target in the concat file
                    concat_fasta_writer.write(line);
                    concat_fasta_writer.write(System.lineSeparator());

                }
                else //one protein read complete, generate decoy and write it
                {

                    //a check validity function is_valid can be used, if we do not want to include certain amino acids e.g, X
                    //however, now it distinguishes among the first '>' marks and the next ones, validity will be turned into true after we get the first
                    // '>' mark. because this is the first annotation, no protein was read yet
                    // so we will start generating decoy from the next time. this time protein is empty.
                    if(validity)//for first annotation this will not execute
                    {
                        //count_AA(protein);// why counting protein again? need to init_composition_hashmap();
                        newline_positions.add(-1);//add a flag

                        check_generate_write(protein, newline_positions, annotation);//the decoy annotation created in the previous iteration
                        decoy_fasta_writer.write(System.lineSeparator());
                        concat_fasta_writer.write(System.lineSeparator());
                    }

                    concat_fasta_writer.write(line);//target annotation
                    concat_fasta_writer.write(System.lineSeparator());

                    protein="";
                    validity=true;


                    /*annotation for the next decoy protein*/

                    //annotation=">"+decoy_tag+line.substring(1, line.length())+System.lineSeparator();

                    //annotation peptideShaker type
                    //int index1=line.indexOf("|");
                    //int index2=line.indexOf("|", index1+1);
                    //annotation=">"+line.substring(1, index1+1)+line.substring(index1+1, index2)+"_"+decoy_tag+line.substring(index2, line.length())+System.lineSeparator();

                    annotation=">"+decoy_tag+line.substring(1, line.length())+System.lineSeparator();

                    //newline positions of the current protein is taken into account while writing the protein,
                    //so we initialize this again, for the new protein to be read next
                    newline_positions=new ArrayList<>();

                }
        }


        if(validity)
        {
            newline_positions.add(-1);//add a flag
            //count_AA(protein);
            check_generate_write(protein, newline_positions, annotation);

        }

        decoy_fasta_writer.write(System.lineSeparator());
        concat_fasta_writer.write(System.lineSeparator());
    }

    /**
     *
     * @param fasta_reader reader object
     * @param generation whether it is generation or statistics
     * @return returns the generated composition HashMap of the fasta file the reader is pointing to
     */
    private HashMap<Character, Double> read_composition(readFASTA fasta_reader, boolean generation)
    {
        de_bruijn_decoy_input = new HashMap<>();//contains kmers? outgoing?

        HashMap<Character, Double> composition = new HashMap<>();
        composition = init_composition_hashmap();

        String line, protein="";
        
        while((line=fasta_reader.readline()) != null)
        {   
            //System.out.println(line);
            if(line.length()>0)
            if(line.charAt(0)!='>')
            {
                protein+=line;//+System.lineSeparator();
        
            }
            else// one protein read complete, skip the annotation and generate kemers and composition for this protein
            {
                composition=count_AA(protein, composition, generation);
                parse_kmer(protein);
                
                protein="";
                
            }
        }
        
        //count and parse the final protein
        composition=count_AA(protein, composition, generation);//update composition with each protein read

        parse_kmer(protein);

        return composition;
    }

    /**
     * Purpose: count number of aminco acids in the protein and updates the composition
     * @param protein the protein that is read
     * @param composition the composition to be updated
     * @param generation if we want to generate new decoy or is this a statistics mode
     * @return the updated composition HashMap
     */
    private HashMap<Character, Double>  count_AA (String protein, HashMap<Character, Double> composition, boolean generation)
    {
        double count_kprp=0.0;
        double count;
        
        char prev_aa=DUMMY;//this variable is needed for statistics, to count number of [KR]P (trypsin) or [WYF]P (chymotrypsin) in data
        
        if(generation)
        {
            if(protein.length()>k)//in generation mode, first kmer's k-1 letters are skipped, because the first k-1 letters are not changed anyways
            {
                prev_aa=protein.charAt(k-1);
                protein=protein.substring(k);
            }
        }

        int length=protein.length();

        for(int i=0; i<length; i++)
        {
            char aa=protein.charAt(i);
            count=composition.get(aa);
            composition.put(aa, count+1);//update composition after each protein read
            
            if(cleavage_site(prev_aa))
            {   
                if(aa=='P')
                {                    
                    count_kprp=count_kprp+1.0;                    
                }
            }
            prev_aa=aa;
        }
        
        count_kprp=count_kprp+composition.get(DUMMY_KP_RP);//to save the target_composition of KP/RP or WP/YP/FP in the protein,
        // needed to know the count of P, in KP/RP or WP/YP/FP
        // we need this to subtract this count from count of 'P',
        composition.put(DUMMY_KP_RP, count_kprp);//if nonspeciifc this should stay 0

        return composition;
    }

    /**
     * Purpose: find kmers present in the target by parsing each protein
     * @param protein the protein to be parsed
     * Effects: the kmers are inserted in the template HashMap de_bruijn_decoy_input including the count of occurrence, which indicates outgoing
     */
    private void parse_kmer(String protein)//create hash for all kmers found in target proteins
    {
        String kmer="";
        int length=protein.length();
        for(int i=0; i<length-(k+1)+1; i++)
        {
            kmer=protein.substring(i, i+(k+1));
        
            if(de_bruijn_decoy_input.containsKey(kmer))
            {
                K_MER count = de_bruijn_decoy_input.get(kmer);
                count.add();  //increases the count
                de_bruijn_decoy_input.put(kmer, count);
            }
            else
                de_bruijn_decoy_input.put(kmer, new K_MER(DUMMY));//initially dummy replacement AA is assigned to all the kmers, outgoing 1
        }
    }

    /**
     * Purpose: counts the total number of Amino Acids in the database by summing composition of all the AAs,
     * needed to compute normalized composition
     * @param composition AA composition of the database
     * @return count of total AAs in the database
     */
    double count_total_AA(HashMap<Character, Double> composition)
    {
        double total=0.0;
        for(Character aa: composition.keySet())
        {
            if(aa!=DUMMY_KP_RP)
            {
                double count=composition.get(aa);
                total+=count;
            }
        }
        
        return total;
    }

    /**
     * Purpose: counts the total number of AAs (without the cleavage AAs (K/R or W/Y/F)) in the database by summing composition of all the AAs,
     * needed to compute normalized composition when cleavage site is kept fixed
     * @param composition AA composition of the database
     * @return count of total AAs (without the cleavage AAs (K/R or W/Y/F)) in the database
     */
    double count_total_AA_wo_cleavage(HashMap<Character, Double> composition)//counts the total number of Amino Acids in the database, needed to compute compostition
    {
        double total=0.0;
        for(Character aa: composition.keySet())
        {
            if(aa!=DUMMY_KP_RP && !cleavage_site(aa))
            {
                double count=composition.get(aa);
                total+=count;
            }
        }
        return total;
    }

    /**
     * Purpose: Debug, print total number of AA
     * @param composition AA composition of the database
     */
    private void print_total_AA(HashMap<Character, Double> composition)
    {        
        double total=count_total_AA(composition);//total number of aa
        //System.out.println("aa count: "+total);
        System.out.printf("AA count: %f\n", total);
        
    }

    /**
     *
     * @param generation if the mode is decoy generation or statistics
     * @return returns the target composition without the fixed AAs (a cleavage site (K/R or W/Y/F), P after [KR](trypsin) or [WYF] (chymotrypsin)
     */
    private HashMap<Character, Double> find_composition_wo_cleavage_AA(boolean generation)
    {
        open_reader_writer(target_fasta_reader, decoy_fasta_writer);

        target_composition=read_composition(target_fasta_reader, generation);
        target_composition_wo_fixedAA=copy_composition_to_hashmap(target_composition, target_composition_wo_fixedAA);//copy target_composition to another variable,
        // the hashmap will contain the original composition
        //generator.calculate_composition();  
        //where is w/o composition? its just copying
        closeReader(target_fasta_reader);        //?rakhbo?
        
        return target_composition_wo_fixedAA;
    }

    /**
     * Purpose: Update the composition HashMap to subtract the counts of P after [KR](trypsin) or [WYF] (chymotrypsin)
     * @param composition AA composition of the database
     * @return the updated composition HashMap
     */
    private HashMap<Character, Double> correct_composition_wo_fixed_AA(HashMap<Character, Double> composition)
    //calculated_composition normalize after subtracting K, R, and KP/RP
    {        
        double total=count_total_AA(composition);//total number of aa

        /////////// to update count of 'P' //////////////
        double count_kprp=composition.get(DUMMY_KP_RP);//KP/RP or WP/YP/FP//
        double count_p=composition.get('P');
        count_p=count_p-count_kprp; //KP/RP or WP/YP/FP is preserved, therefore do not need this count in target_composition
        
        if(count_p <0) {
            System.exit((int)count_p);
        }
        
        composition.put('P', count_p); //updated 'P' excluding the number of KP/RPs
        
        return composition;
    }
    
    
    
    //ekta sub korei abar shob new calculate kortesi, dorkar nai. update kora jabe hishab korlei?

    /**
     * purpose: the normalized composition without considering the cleavage amino acids of corresponding enzyme ([KR](trypsin) or [WYF] (chymotrypsin))
     * @param composition AA composition of the database
     * @return returns the normalized HashMap<Character, Double> type, which contains the normalized composition
     */
    private HashMap<Character, Double> calculate_normalized_composition_wo_cleavage_AA(HashMap<Character, Double> composition) //calculated_composition normalize after subtracting K, R, and KP/RP
    {        
        
        double total=count_total_AA_wo_cleavage(composition);//total number of aa
        HashMap<Character, Double> normalized = normalize_composition(total, composition);//for adaptive normalize pore korle valo hobe, alada ekta variable a, tahole each time subtract and normalize
        
        return normalized;
    }

    /**
     * Purpose: return the normalized composition, given the composition of AAs in the database
     * @param total total number of AAs
     * @param composition  AA composition of the database
     * @return the normalized composition according to the provided total
     */
    private HashMap<Character, Double> normalize_composition(double total, HashMap<Character, Double> composition)
    {        
        HashMap<Character, Double> normalized = new HashMap<>();
        
        for (Character aa : composition.keySet()) {
            double count=composition.get(aa); //count is a whole number
            double frequency=total==0.0? 0.0:count/total;   //divide by total number of aa, to get the target_composition/frequency of each aa
            normalized.put(aa, frequency); //put frequency which is in the range of 0 to 1
        }
        
        return normalized;
    }

    /**
     * Purpose: to calculate the cumulative probability, so that random AA can be selected using the probability distribution
     * @param composition AA composition of the database
     * @return ArrayList<AminoAcid>  containing each AA and their cumulative probability
     */
    private ArrayList<AminoAcid> calculate_cum_probability(HashMap<Character, Double> composition)//calculate cumulative, sorted in increasing order
    {  
        ArrayList<AminoAcid> cumulative_prob =new ArrayList<>();//new calculation start, check first time calculation too?
        cumulative_prob = copy_composition(cumulative_prob, composition);//copy of probability of aa, Amino Acid needed for lenght and mass distribution, do not need this object anymore
        cumulative_prob = sort_composition(cumulative_prob);// sort min to max
        
        AminoAcid prevAA=new AminoAcid();

        double prob=prevAA.probability;
        
        for(AminoAcid AA: cumulative_prob)
        {
            prob=prevAA.probability;
            double comp=AA.composition;//current aa composition
            if(comp!=0 && ((AA.aa!='X'&& AA.aa!='K'&& AA.aa!='R' && ENZYME==TRYPSIN)||(AA.aa!='X'&& AA.aa!='Y' && AA.aa!='W' && AA.aa!='F' && ENZYME==CHYMOTRYPSIN)) || (ENZYME==NONSPECIFIC))//do not need to check K/R as not copied
            {
                prob=comp+prevAA.probability;
                prevAA=AA;
            }
            
            AA.probability=prob;//probability updated, composition fixed
            
        }

        return cumulative_prob;
    }

    /**
     * Purpose: Copy source src_composition (HashMap) to destination dst_cp_probability_AA (ArrayList of AminoAcid) except the fixed AAs
     * fixed AA: a cleavage site (K/R or W/Y/F), P after [KR](trypsin) or [WYF] (chymotrypsin)
     * @param dst_cp_probability_AA destination to where the composition will be copied
     * @param src_composition source composition from where composition will be copied
     * @return destination dst_cp_probability_AA ArrayList of AminoAcid is returned
     */
    ArrayList<AminoAcid> copy_composition(ArrayList<AminoAcid> dst_cp_probability_AA, HashMap<Character, Double> src_composition)//copy src_composition to a list of AminoAcids
    {        
        for(Character aa: src_composition.keySet())
        {
            if(aa!=DUMMY_KP_RP && !cleavage_site(aa))//every letter in cumulative_sum_prob_AA except for the fixed ones (a cleavage site (K/R or W/Y/F), P after [KR](trypsin) or [WYF] (chymotrypsin)
            {
                double comp= src_composition.get(aa);
                dst_cp_probability_AA.add(new AminoAcid(aa, comp));
            }
        }
        return dst_cp_probability_AA;
    }


    /**
     * Purpose: sort the AminoAcid ArralyList according to the AAs composition
     * @param sprobability_AA the ArraList to be sorted
     * @return
     */
    ArrayList<AminoAcid> sort_composition(ArrayList<AminoAcid> sprobability_AA)
    {
        sprobability_AA.sort(AminoAcidComparator);
        return sprobability_AA;
    }

    /**
     * comparator to sort ArrayList, ascending order
     */
    public static Comparator<AminoAcid> AminoAcidComparator = (AminoAcid a1, AminoAcid a2) -> {
        //ascending order
        if (a1.composition>a2.composition)
            return -1;
        else if (a1.composition<a2.composition)
            return 1;
        else
            return 0;
    };


    /**
     * Purpose: print for debugging purpose, without the fixed AAs
     * fixed AA: a cleavage site (K/R or W/Y/F), P after [KR](trypsin) or [WYF] (chymotrypsin)
     * @param composition the AA composition to be printed
     */
    void print_composition_wo_clevage_AA(HashMap<Character, Double> composition)
    {
        String log = "";
        double total = 0.0;
        for(Character aa: composition.keySet())
        {
            double count=composition.get(aa);

            if(aa=='_') {
                //System.out.println("KP/RP: " + count);
                log = log + "KP/RP: " + count + "\n";
            }
            else
            {
                if(!((ENZYME==TRYPSIN&&(aa=='K'||aa=='R'))||(ENZYME==CHYMOTRYPSIN&&(aa=='W'||aa=='Y'||aa=='F'))))
                    total += count;
                //System.out.println(aa+": "+count);
                log = log+aa+": "+count+"\n";
            }
        }
        log = log+"\nTotal composition (W/O KP/RP or WP/YP/FP): "+total;
        //System.out.println("Total composition (W/O KP/RP or WP/YP/FP): "+total);
    }

    /**
     * Purpose: print for debugging purpose, without the fixed AAs
     * @param composition composition the AA composition to be printed
     */
    void print_composition(HashMap<Character, Double> composition)
    {
        String log = "";
        double total = 0.0;
        for(Character aa: composition.keySet())
        {
            double count=composition.get(aa);
            
            
            if(aa=='_')
            {
                //System.out.println("KP/RP: "+count);
                log = log+"KP/RP: "+count+"\n";
            }
            
            else
            {
                total += count;
                //System.out.println(aa+": "+count);
                log = log+aa+": "+count+"\n";
            }
        }
        
        //System.out.println("Total composition (W/O KP/RP or WP/YP/FP): "+total);
        log = log+"\nTotal composition (W/O KP/RP or WP/YP/FP): "+total;
        
        LOGGER.finer(log);
    }

    /**
     * Purpose: print probability, for debugging purpose
     * @param cumsum_prob the ArrayList of containing the AA probability
     */
    void print_probability(ArrayList<AminoAcid> cumsum_prob) //for debugging and testing purpose
    {
        String log = "Cumulative probability\n";

        double last_AA_prob = 0;
        for(AminoAcid AA: cumsum_prob)
        {
            if(AA.probability!=0)
                last_AA_prob = AA.probability;
            log = log+AA.aa+" "+AA.probability+"\n";
        }
        

        log = log+" \nlast probability = "+last_AA_prob+"\n";
        
        LOGGER.finer(log);
    }

    /**
     * Purpose: write decoy annotation, check the target protein, generate the decoy
     * @param protein target protein
     * @param newline_positions position of the newlines in the target, needed only for formatting, supplied to the method "generate_decoy"
     * @param annotation generated decoy annotation
     */
    void check_generate_write(String protein, ArrayList <Integer> newline_positions, String annotation)
    {
        if(!protein.equals(""))
        {
            decoy_fasta_writer.write(annotation);//annotation for the protein to be read next
            concat_fasta_writer.write(annotation);//annotation for the protein to be read next

            for(int i=0; i<k; i++)//keep the first k same
            {
                decoy_fasta_writer.write(protein.charAt(i));
                concat_fasta_writer.write(protein.charAt(i));
            }

            generate_decoy(protein, newline_positions);
            //write(System.lineSeparator());
        }
    }


    /**
     * Purpose: write deBruijn decoy template to a file
     * @param de_bruijn deBruijn decoy template
     */
    void write_decoy_template(HashMap<String, K_MER> de_bruijn)
    {
        int i=1;
        for( String kmer: de_bruijn.keySet())
        {
            de_bruijn_writer.write(kmer+" -> "+de_bruijn.get(kmer).decoy+" ");
            if(i%maximum_entry_per_line==0)
                de_bruijn_writer.write(System.lineSeparator());
            i++;
        }
        System.out.println();
    }


    /**
     * Purpose: select the deBruijn template to be used
     * @param use_input if true, use the deBruijn template generated considering target AA composition,
     *                  if false, use the deBruijn template generated usimng random AA
     */
    void select_de_bruijn(boolean use_input)
    {
        if(use_input)
            de_bruijn=de_bruijn_decoy_input;
        else
            de_bruijn=de_bruijn_decoy_random;
    }
    
    void getLoggerFile()
    {
        
        Handler consoleHandler = null;
        Handler fileHandler  = null;
        try{
            //Creating consoleHandler and fileHandler
            consoleHandler = new ConsoleHandler();
            fileHandler  = new FileHandler("./deBruijn.log");
             
            //Assigning handlers to LOGGER object
            LOGGER.addHandler(consoleHandler);
            LOGGER.addHandler(fileHandler);
             
            //Setting levels to handlers and LOGGER
            consoleHandler.setLevel(Level.ALL);
            fileHandler.setLevel(Level.ALL);
            LOGGER.setLevel(Level.ALL);
             
            //LOGGER.config("Configuration done.");
             
            //Console handler removed
            LOGGER.removeHandler(consoleHandler);
             
            LOGGER.log(Level.FINE, "Finer logged");
        }catch(IOException exception){
            LOGGER.log(Level.SEVERE, "Error occur in FileHandler.", exception);
        }
         
        LOGGER.finer("Finest example on LOGGER handler completed.");
    }


    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        long startTime = System.currentTimeMillis();
        boolean use_input=true;
               
        //String inputfilenName="uniprot-proteomeUP000005640_20180206.fasta";//"";//"uniprot-mus-musculus.fasta";//"uniprot-human-reviewed-Nov-2017-Trypsin.fasta";//"swiss_prot_all_uncompressed_08082018.fasta";//"ups.fasta";//


        ENZYME = TRYPSIN;//NONSPECIFIC;//


        Args_Parser parser = new Args_Parser();
        parser.parse_args(args);
        parser.check_args();

        String inputfilenName = parser.get_inputfilename();
        String inputDirectory = parser.get_inputDirectory();
        //String inputfilename = parser.get_inputfilename();
        deBruijn.k = parser.get_k();
        ENZYME = parser.get_ENZYME();

        String prefix="e"+ENZYME+"k"+k+"DBD-";
        
        //give the option to select output file name
        String outputfilename=parser.get_outputfilename();
        deBruijn generator = new deBruijn(inputDirectory, inputfilenName, outputfilename, k);//DEBRUIJN_DECOY//
        

        //Step 1 CALCULATE COMPOSITION
        
        /////////////////////to calculate the target target_composition wo kr kp rp (as they will be fixed)//////////////////////////

        generator.target_composition_wo_fixedAA=generator.find_composition_wo_cleavage_AA(true);//to generate the decoy
        
        generator.normalized_composition=generator.calculate_normalized_composition_wo_cleavage_AA(generator.target_composition_wo_fixedAA);//calculated_composition normalize after subtracting K, R, and KP/RP
        
        
        //System.out.println("TARGET COMPOSITION (W/O X and K, R, KP, RP or W, Y, F, WP, YP, RP)");
        //generator.print_composition_wo_clevage_AA(generator.normalized_composition);    //no need to print this, only for debugging purpose    //more than one because of K/R, W/Y/F (fixed locations)
        
        /////////////////////////////////////////////
               
        
       //Step2 GENERATE DEBRUIJN DECOY
        generator.adaptive_composition=generator.copy_composition_to_hashmap(generator.target_composition_wo_fixedAA, generator.adaptive_composition);//adaptive use korle tai copy kore rakhlam
        generator.adaptive_composition=generator.correct_composition_wo_fixed_AA(generator.adaptive_composition);
        generator.normalized_composition=generator.calculate_normalized_composition_wo_cleavage_AA(generator.adaptive_composition);
        cummulative_sum_prob_AA =generator.calculate_cum_probability(generator.normalized_composition);//cummulitive probabilty, to do weighted sum sampling

        //generator.calculate_cum_probability(generator.target_composition_wo_fixedAA);//cummulitive probabilty, to do weighted sum sampling
        generator.print_probability(cummulative_sum_prob_AA); //print cumsum probability, for debugging purpose //without normalization?


        // random generation
        generator.generate_decoy_template_random(generator.de_bruijn_decoy_random);
        //generator.print_decoy_template(generator.de_bruijn_decoy_random);

        // generation considering input composition
        generator.adaptive_composition=generator.generate_decoy_graph_input_random(generator.de_bruijn_decoy_input, generator.adaptive_composition, cummulative_sum_prob_AA);
        //generator.print_decoy_template(generator.de_bruijn_decoy_input);
        generator.normalized_composition=generator.calculate_normalized_composition_wo_cleavage_AA(generator.adaptive_composition);
        generator.print_composition_wo_clevage_AA(generator.normalized_composition);


        generator.select_de_bruijn(use_input);

        generator.open_reader_writer(generator.target_fasta_reader, generator.decoy_fasta_writer);
        generator.read_generate_write();

        //generator.print_decoy_template(generator.de_bruijn);//print the selected one
        generator.write_decoy_template(generator.de_bruijn);//write the selected one


        //////////////////////////////////
       
        generator.closeFile(generator.target_fasta_reader, generator.decoy_fasta_writer);
        
        //////////////////////////////////
  
        long endTime = System.currentTimeMillis();
        
        long timeElapsed = endTime - startTime;

        
        System.out.println("DeBruijn Decoy Generated, total time: "+timeElapsed+" ms");
        
    }
    
}
