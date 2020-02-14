/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deBruijn;

/**
 *
 * @author jmmoosa
 */

import java.io.*;

import java.util.HashMap;

public class readFASTA {
    
    public HashMap<String, Peptide> digested_u_peptides;
    
    
    String inputfilename;
    
    FileReader fileReader;
    BufferedReader bufferedReader;
        
    readFASTA(String inputfilename, HashMap<String, Peptide> digested_u_peptides)
    {
        
        this.inputfilename=inputfilename;
        
        this.digested_u_peptides=digested_u_peptides;
        //openReader();
        //openWriter();
    }
    
    readFASTA(String inputfilename)
    {
        this.inputfilename=inputfilename;
        
    }
    
    void openReader()
    {
        try {
            
            fileReader = new FileReader(inputfilename);
            bufferedReader = new BufferedReader(fileReader);
            
        }
        catch(FileNotFoundException ex) {
            System.out.println(
                "Unable to open file '" + 
                inputfilename + "' to parse");  
            System.exit(1);
        }
        
    }
    
    String readline()
    {
        String line;
        boolean validity=true;
        try
        {
            
            if((line=bufferedReader.readLine()) != null) {
            
                return line;
            }
            
        }
        catch(IOException ex) {
            System.out.println(
                "Error reading file '" 
                + inputfilename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }

       return null;
    }
    
    
   
    void closeReader()
    {
        try
        {            
            bufferedReader.close(); 
            fileReader.close();

        }
        catch(IOException ex) {
            System.out.println(
                "Error reading file '" 
                + inputfilename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }
    }
    
    
}
