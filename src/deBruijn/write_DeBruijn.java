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

public class write_DeBruijn {
    
    String outputfilename;
    
    FileWriter fileWriter;
    BufferedWriter writer;
        
    write_DeBruijn(String outputfilename)
    {
        this.outputfilename=outputfilename;
        
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
    
    
    void write (String line)
    {
        try {
            
            writer.append(line);
        }
        catch(IOException ex) {
            System.out.println(ex+
                "Error writing decoy deBruijn template in file '" 
                + outputfilename + "'");                  
            // Or we could just do this: 
            // ex.printStackTrace();
        }
        
    }
    
    void write (char decoy)
    {
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


