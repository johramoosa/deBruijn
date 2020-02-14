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
public class AminoAcid {
    char aa;
    double composition, probability;
    
    String single_letter_code, three_letter_code;
    double residue_mass, average_mass, immonium_mass;//charge for immo?
    
    AminoAcid(String single_letter_code, String three_letter_code, double residue_mass, double average_mass, double immonium_mass)
    {
        this.residue_mass=residue_mass;
        this.single_letter_code=single_letter_code;
        this.three_letter_code=three_letter_code;
        this.average_mass=average_mass;
        this.immonium_mass=immonium_mass;
    }
    
    AminoAcid(String single_letter_code, String three_letter_code, double residue_mass, double average_mass)
    {
        this.residue_mass=residue_mass;
        this.single_letter_code=single_letter_code;
        this.three_letter_code=three_letter_code;
        this.average_mass=average_mass;
        this.immonium_mass=0.0;
    }
    
    
    AminoAcid()
    {
        composition=0.0;
        probability=0.0;
    }
    AminoAcid(char aa, double composition)
    {
        this.aa=aa;
        this.composition=composition;
        
        probability=0.0;
    }
    
}
