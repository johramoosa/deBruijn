package deBruijn;

public class Option {
    String key;
    int min_len=1;//min length of the argument
    String type;//the type of the value
    private boolean optional;//if this option is optional
    private boolean value_required;//if providing a value is required if this option is given
    private boolean provided_option;//if this option is provided
    private boolean provided_value;//if any value is provided for this option
    boolean skip;
    String value;

    Option()
    {
        min_len = 1;
        type = "int";
    }
    Option(int min_len, String type)
    {
        this.min_len = min_len;
        this.type = type;
        provided_option = false;
        provided_value = false;
    }

    Option(String key, int min_len, String type)
    {
        this.key = key;
        this.min_len = min_len;
        this.type = type;
        optional = true;
        provided_option = false;
        provided_value = false;
        value_required = false;
    }

    Option(String key, int min_len, String type, boolean value_required, boolean optional)
    {
        this.key = key;
        this.min_len = min_len;
        this.type = type;
        this.optional = optional;
        provided_option = false;
        provided_value = false;
        this.value_required = value_required;
    }

    Option(String key, int min_len, String type, boolean value_required)
    {
        this.key = key;
        this.min_len = min_len;
        this.type = type;
        optional = true;
        provided_option = false;
        provided_value = false;
        this.value_required = value_required;
    }

    void setProvidedOption()
    {
        provided_option = true;
    }

    boolean get_provided_value()
    {
        return provided_value;
    }
    void setProvidedValue()
    {
        if(provided_option)
            provided_value = true;
    }

    void print()
    {
        if(key.equals(Args_Parser.HELP_key))
        {
            System.out.println("--input <*.fasta> \n" +
                    "[--output <*.fasta>]\n" +
                    "[--k <integer value greater than 0>] (default: 2)\n" +
                    "[--enzyme 0/1/2] (0=non specific, 1=trypsin (default), 2=chymotrypsin)\n" +
                    "[--compositoin 0/1/2] (0=not adaptive, target composition not considered, 1=not adaptive, target composition considered, 2=adaptive, target composition consodered (default))");
        }
    }

    void set_value(String value)
    {
        this.value= value;
        setProvidedValue();

    }

    boolean check_valid()
    {
        boolean is_Valid = (provided_option && value_required && provided_value && check_value())||
                (provided_option&&((value_required&&provided_value&&check_value())||(!value_required&&!provided_value)));
        return is_Valid;
    }

    boolean check_required()//return false if required input not given
    {
        boolean required = (!optional&&provided_option)||optional;
        return required;
    }

    boolean check_skip ( )
    {
        skip = (optional && !provided_option);
        return skip;
    }

    void print_error()
    {
        if(!optional&&!provided_option)
            System.out.println("Required Option "+key+" not provide");
        else if(optional&&provided_option&&value_required&&!provided_value)
            System.out.println("Option "+key+" requires a value");
        else if(optional&&provided_option&&value_required&&provided_value&&!check_value()) {
            System.out.print("Wrong value. Option " + key + " requires a value of the format: ");
            if(key.equals(Args_Parser.INPUT_key)||key.equals(Args_Parser.OUTPUT_key))
                System.out.println("<filename>.fasta");
            else if(key.equals(Args_Parser.COMPOSITION_key)||key.equals(Args_Parser.ENZYME_key))
            {
                System.out.println("Integer: 0/1/2");
            }
            else if(key.equals(Args_Parser.K_key))
            {
                System.out.println("Positive Integer");
            }
        }


    }
    boolean check_value()
    {
        if(key.equals(Args_Parser.INPUT_key))
        {
            if(value.endsWith(".fasta") && value.length()>min_len)
            {
                return true;
            }
        }
        else if(key.equals(Args_Parser.HELP_key))
        {
            return true;
        }
        else if(key.equals(Args_Parser.OUTPUT_key))
        {
            if(value.endsWith(".fasta") && value.length()>min_len)
            {
                return true;
            }
        }
        else if(key.equals(Args_Parser.COMPOSITION_key))
        {
            try{
                int v = Integer.parseInt(value);
                if(v>0 && v<3)
                    return true;
            }catch(NumberFormatException e){
                return false;
            }
        }
        else if(key.equals(Args_Parser.K_key))
        {
            try{
                int v = Integer.parseInt(value);
                if(v>0)
                    return true;
            }catch(NumberFormatException e){
                return false;
            }
        }
        else if(key.equals(Args_Parser.ENZYME_key))
        {
            try{
                int v = Integer.parseInt(value);
                if(v>0 && v<3)
                    return true;
            }catch(NumberFormatException e){
                return false;
            }
        }


        return false;
    }
    boolean check_type(Object o)
    {
        if(o instanceof Integer)
        {
            if(type.equals("int")) {
                return true;
            }
        }
        else if(o instanceof String)
        {
            if(type.equals("String")) {
                return true;
            }
        }
        else if(o instanceof Boolean)
        {
            if(type.equals("boolean")) {
                return true;
            }
        }

        return false;
    }
}
