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
//to write the decoy sequence
import java.io.*;
import java.io.FileWriter;
import java.util.HashMap;

public class writeFASTA {
    
    public HashMap<String, Peptide> digested_u_peptides;
    
    
    String outputfilename;
    
    FileWriter fileWriter;
    BufferedWriter writer;
        
    writeFASTA(String inputfilename, HashMap<String, Peptide> digested_u_peptides)
    {
        
        this.digested_u_peptides=digested_u_peptides;
        //openReader();
        //openWriter();
    }
    
    writeFASTA(String outputfilename)
    {
        
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
