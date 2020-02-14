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
import java.io.FileWriter;
import java.util.HashMap;

public class read_writeFASTA {
    
    public HashMap<String, Peptide> digested_u_peptides;
    
    
    String inputfilename, outputfilename;
    
    FileReader fileReader;
    BufferedReader bufferedReader;
    FileWriter fileWriter;
    BufferedWriter writer;
        
    read_writeFASTA(String inputfilename, HashMap<String, Peptide> digested_u_peptides)
    {
        this.inputfilename=inputfilename;
        
        this.digested_u_peptides=digested_u_peptides;
        //openReader();
        //openWriter();
    }
    
    read_writeFASTA(String inputfilename, String outputfilename)
    {
        this.inputfilename=inputfilename;
        this.outputfilename=outputfilename;
        //this.digested_u_peptides=digested_u_peptides;
        //openReader();
        //openWriter();
    }
    
    void openWriter()
    {
        try {
            
            fileWriter=new FileWriter(outputfilename);
            writer = new BufferedWriter(fileWriter);
            
        }
        catch(IOException ex) {
            System.out.println(
                "Unable to open file '" + 
                outputfilename + "' to write random sequence");                
        }        
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
    
    
    void write (String protein)
    {
        //System.out.println("String write..................."+protein);
        try {
            
            writer.append(protein);
        }
        catch(IOException ex) {
            System.out.println(ex+
                "Error writing decoy sequence string "+protein+" in file '" 
                + outputfilename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }
        
    }
    
    void write (char decoy)
    {
        //System.out.println("char write..................."+decoy);
        try {
            
            writer.append(decoy);
        }
        catch(IOException ex) {
            System.out.println(ex+
                "Error writing decoy AA "+decoy+" in file '" 
                + outputfilename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }
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
    
    void closeWriter()
    {
        try {
            writer.flush();
            writer.close();
            fileWriter.close();
        }
        catch(IOException ex) {
            System.out.println(ex+
                " Error closing file '" 
                + outputfilename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }
    }
    
    
    
    
        
    
}
