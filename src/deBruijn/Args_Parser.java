package deBruijn;

import java.util.HashMap;

public class Args_Parser {


    private HashMap<String, String> hash_key_hash;

    private HashMap<String, Option> key_hash;

    static final String HELP_key = "help";
    static final String INPUT_key = "input";
    static final String OUTPUT_key = "output";
    static final String K_key = "k";
    static final String ENZYME_key = "enzyme";
    static final String COMPOSITION_key = "composition";


    private String inputfilename="", outputfilename="";
    private Integer ENZYME=1;
    private Integer k=2;
    private Boolean consider_target_composition=true, isAdaptive=true;


    Args_Parser()
    {
        hash_key_hash = new HashMap<>();
        init_Hash_Hash_Keys();
        init_Hash_Keys();
    }

    /**
     * Purpose: to initialize options value pairs
     */
    private void init_Hash_Keys()
    {
        key_hash = new HashMap<>();

        Option opt = new Option(HELP_key, 0, "null");//optional, no value required
        key_hash.put(HELP_key, opt);

        opt = new Option(INPUT_key, 6, "String", true, false);//only one that is a must required option
        key_hash.put(INPUT_key, opt);

        opt = new Option(OUTPUT_key, 1, "String", true);//optional
        key_hash.put(OUTPUT_key, opt);

        opt = new Option(K_key, 1, "int", true);//optional
        key_hash.put(K_key, opt);

        opt = new Option(ENZYME_key, 1, "int", true);//optional
        key_hash.put(ENZYME_key, opt);

        opt = new Option(COMPOSITION_key, 1, "int", true);//or 0/1/2 0 = random, no composition
        // 1 = composition non adaptive, 2 = adaptive composition (default) combine this two
        key_hash.put(COMPOSITION_key, opt);
    }

    /**
     * Name: private void init_Hash_Hash_Keys()
     * Input: None
     * Output: None
     * Effect: initializes the hash_key_hash variable
     * Purpose: To initialize the Hash key mapping Hash, so that multiple keys can point to same value, i.e., so that multiple keywords or keys can be used for one single option
     */
    private void init_Hash_Hash_Keys()
    {
        hash_key_hash = new HashMap<>();
        //HELP

        hash_key_hash.put("help", HELP_key);
        hash_key_hash.put("h", HELP_key);

        //INPUT, filename.fasta (.fasta na thakle incorrect file type)//need input key before output key

        hash_key_hash.put("input", INPUT_key);
        hash_key_hash.put("in", INPUT_key);
        hash_key_hash.put("i", INPUT_key);

        //OUTPUT, filename.fasta (.fasta na thakle incorrect file type)
        hash_key_hash.put("output", OUTPUT_key);
        hash_key_hash.put("out", OUTPUT_key);
        hash_key_hash.put("o", OUTPUT_key);

        //K, integer number the higher the slower, min:1 , max:5

        hash_key_hash.put("k", K_key);//use - only?

        //ENZYME, nonspecific = 0, trypsin = 1, chymotrypsin = 2

        hash_key_hash.put("enzyme", ENZYME_key);
        hash_key_hash.put("e", ENZYME_key);

        //USE COMPOSITION, boolean

        hash_key_hash.put("composition", COMPOSITION_key);
        hash_key_hash.put("c", COMPOSITION_key);

        //USE ADAPTIVE, boolean

        //hash_key_hash.put("adaptive", ADAPTIVE_key);
        //hash_key_hash.put("a", ADAPTIVE_key);
    }

    void check_args()
    {
        for(String key: key_hash.keySet())
        {
            Option opt = key_hash.get(key);
            //System.out.print(key+"...key");
            if(opt.check_required()) {
                if(opt.check_skip())
                {
                    //this option skipped
                }
                else {
                    if (opt.check_valid( )) {
                        if (key.equals(INPUT_key)) {
                            inputfilename = opt.value;
                            System.out.println("Input filename: " + inputfilename);
                        }
                        else if (key.equals(OUTPUT_key)) {
                            outputfilename = opt.value;
                            System.out.println("Output filename: " + outputfilename);
                        }
                        else if (key.equals(HELP_key)) {
                            opt.print( );
                            System.exit(0);
                        }
                        else if (key.equals(ENZYME_key)) {
                            ENZYME = Integer.parseInt(opt.value);
                            System.out.println("ENZYME: " + ENZYME);
                        }
                        else if (key.equals(K_key)) {
                            k = Integer.parseInt(opt.value);
                            System.out.println("k: " + k);
                        }
                        else if (key.equals(COMPOSITION_key)) {
                            int compo = Integer.parseInt(opt.value);
                            if (compo == 0) {
                                consider_target_composition = false;
                                isAdaptive = false;
                            }
                            else if (compo == 1) {
                                consider_target_composition = true;
                                isAdaptive = false;
                            }
                            else if (compo == 2)//default
                            {
                                consider_target_composition = true;
                                isAdaptive = true;

                            }

                            System.out.println("consider_target_composition: " + consider_target_composition);
                            System.out.println("isAdaptive: " + isAdaptive);
                        }
                    }
                    else {
                        opt.print_error( );
                        System.exit(1);
                    }
                }
            }
            else
            {
                //System.out.println("requirement error?");
                opt.print_error();
                System.exit(1);
            }
        }


        Option opt = key_hash.get(Args_Parser.OUTPUT_key);
        String prefix="e"+ENZYME+"k"+k+"DBD-";
        if (opt.check_skip()) {
            outputfilename = prefix+inputfilename;
            //System.out.println(outputfilename);
        }

    }
    void parse_args(String[] args)
    {
        int number_of_options = 7;


        String option="h", value, key="";
        Option opt = new Option();
        boolean error_input_format = false;
        boolean error_missing_option = false;
        boolean error_unknown_option = false;
        boolean error_missing_value = false;
        boolean expect_option = true;
        boolean expect_value = false;
        //args gula dekho
        int number_args = args.length;
        for(int i = 0; i<number_args; i++)
        {
            String str = args[i].toLowerCase();

            if(expect_option) {
                if (str.startsWith("-")) {
                    //generate the key to the key for options hash
                    if (str.startsWith("--"))//double hyphen, expecting a word
                    {
                        option = str.substring(2);
                    }
                    else//single hyphen, expecting only a letter after that
                    {
                        if (str.length( ) > 2)
                            error_input_format = true;
                        else if (str.length( ) < 2)
                            error_missing_option = true;
                        else if (str.length( ) == 2) {
                            key = str.substring(1);
                            option = hash_key_hash.get(key);
                        }
                    }

                    //use the key to get the attributes of that option
                    if (hash_key_hash.containsKey(option)) {
                        key = hash_key_hash.get(option);
                        opt = key_hash.get(key);
                        opt.setProvidedOption();

                        if(!key.equals(HELP_key)) {
                            expect_value = true;
                            expect_option = false;
                        }
                        else//help option given
                        {
                            opt.setProvidedOption();
                            key_hash.put(key, opt);
                            expect_value = false;
                            expect_option = false;//help message will be printed. other options and values ignored
                        }
                    }
                    else{
                        error_unknown_option = true;
                        System.exit(1);
                    }
                }
                else
                {
                    error_input_format = true;
                    System.out.println("Error in input: Expected option found value");
                    System.exit(1);
                }
            }
            else if(expect_value)
            {
                if (str.startsWith("-")) {
                    error_missing_value = true;
                    System.out.println("Error: missing value for the option: "+key);
                    System.exit(1);
                }

                value = str;
                opt.set_value(value);//set in the hash?

                expect_value = false;
                expect_option = true;
            }
        }
    }

    String get_inputfilename()
    {
        return inputfilename;
    }

    Integer get_k()
    {
        return k;
    }

    Integer get_ENZYME()
    {
        return ENZYME;
    }

    String get_outputfilename()
    {
        return outputfilename;
    }



}
