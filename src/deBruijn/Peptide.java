/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deBruijn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Random;

/**
 *
 * @author jmmoosa
 */
public class Peptide {
    
    String pepid;
    long id;
    String sequence, shuffled;
    double quality_score;
    int length;
    
    double pepmass;
    String protein_annotation;
    int protein_id;
    int start_index;
    
    int count;
    //massValue ION_mass[];
    
    HashMap<Integer, ArrayList<Integer>> digested_peptides_location;
    
    Peptide()
    {
        sequence="";
        shuffled="";
        quality_score=Double.NEGATIVE_INFINITY;
        pepmass=0;
        //protein_id=-1;
        //count=1;
    }
    
    void initPeptide(long id, String sequence)
    {
        this.digested_peptides_location=new HashMap<>();
        //ION_mass=new massValue[sequence.length()];
        count=1;
        this.id=id;
        this.sequence=sequence;
        quality_score=0.0;
        length=sequence.length();
        calculate_pepmass();
        //calculate_ionmass();
    }
    
    void initPeptide(long id, String sequence, double pepmass)
    {
        initPeptide(id, sequence);
        this.pepmass=pepmass;
    }
    
    Peptide(long id, String sequence)
    {
        initPeptide(id, sequence);
    }
    
    Peptide(long id, String sequence, int protein_id, int start_index)
    {
        initPeptide(id, sequence);
        this.protein_id=protein_id;
        this.start_index=start_index;        
    }
    
    Peptide(long id, String sequence, double pepmass)
    {
        initPeptide(id, sequence, pepmass);
    }
    
    Peptide(long id, String sequence, String protein_annotation)
    {
        initPeptide(id, sequence);
        this.protein_annotation=protein_annotation;
    }
    
    Peptide(int id, String sequence, int protein_id)
    {
        initPeptide(id, sequence);
        this.protein_id=protein_id;
    }
    
    void print()
    {
        //System.out.println(id+" "+quality_score+" "+sequence);
        System.out.println(id+": "+count+" "+pepmass+" "+sequence+", "+digested_peptides_location.values());//+" "+CS682A4Johra.protein_annotation.get(protein_id).substring(0, 10)+" id:"+protein_id);
    }
    
    
    int check_pattern(char aa, char next_aa)//returns true if pattern is not k/r, if kp/rp return true, so comp_stat will find mismacth in number of k and r
    {
        //assert kmer prev
        int skip_count=0;
        
        /*if(next_aa=='R' || next_aa=='K')
        {
            count++;
        }               
        else */
        if(aa=='R' || aa=='K')
        {
            skip_count=1;
            if(next_aa=='P')
                skip_count++;
        }
        return skip_count;
    }
    
    int check_pattern(char aa)//returns true if pattern is not k/r, if kp/rp return true, so comp_stat will find mismacth in number of k and r
    {
        //assert kmer prev
        int skip_count=0;
        
        if(aa=='R' || aa=='K')
        {
            skip_count=1;
            
        }
        return skip_count;
    }
    
    
    void shuffle()
    {
        boolean last_kp_rp=false;
        Random rand;
        rand = new Random();
        shuffled=sequence;
        //need to do this
        //last a K/R thakle dont change
        //dont change any K/R position, and P if age K/R thake
        
        //System.out.println("................"+sequence);
        int length=sequence.length();
        int count_KRP=0;
        
        ArrayList <Character> characters= new ArrayList<>();
        TreeMap<Integer, Character> fixed_positions=new TreeMap<>();
        
        StringBuilder shuffled_c=new StringBuilder (length);
        
        char curr_AA, next_AA, prev_AA;
        curr_AA=sequence.charAt(0);
        next_AA=curr_AA;
        prev_AA=curr_AA;
        for(int i=1; i<length; i++)
        {
            next_AA=sequence.charAt(i);
            //System.out.println(i+"-"+curr_AA+".."+next_AA);
            prev_AA=curr_AA;
            int c;
            if((c=check_pattern(curr_AA, next_AA))>0)
            {
                //System.out.println(".........%%%%%%......."+c);
                count_KRP+=c;
                fixed_positions.put(i-1, curr_AA);
                if(c>1)
                    fixed_positions.put(i, next_AA);//P add hoye gese ekhane
                
                //set curr_AA here, and i
                i++;//RP/KP paisi
                if(c>1 && i<length)
                {
                    curr_AA=sequence.charAt(i);
                    prev_AA=curr_AA;
                    //i++;
                }
                else
                    curr_AA=next_AA;
                
            }
            else
            {
                characters.add(curr_AA);  
                prev_AA=curr_AA;
                curr_AA=next_AA;
            }
        }
        
        //System.out.println(prev_AA+"....********......"+curr_AA);
        
        if(check_pattern(curr_AA)==1)//for the last K/R
        {
            //System.out.println("-"+curr_AA+".."+next_AA);
            fixed_positions.put(length-1, curr_AA);
            count_KRP++;
        }
        else if(check_pattern(prev_AA, curr_AA)>1)//need to check if the last is KP/RP
        {
            //do nothing
            
        }
        else
            characters.add(curr_AA); 
        //System.out.println("...."+curr_AA+"********......"+next_AA);
        int count_c=length-count_KRP;//-1 because i want index, -1 er karone last 2ta NON_SPECIFIC_AA same thaktesilo, 
        
        //int i=0;
        int j=0;
        //int nl=newline_positions.get(j);
        //System.out.println(characters.size()+"................"+count_c);
        //System.out.println(sequence+"-"+characters+" : "+fixed_positions.toString());
        
        while(!characters.isEmpty())
        {   
            //first NON_SPECIFIC_AA can not be P ensure hereee?
           //int random_index= rand.nextInt(characters.size());
            int random_index= 0;
            if (characters.size()>1)
                random_index=rand.nextInt(count_c);
            count_c--;
            
                        
            shuffled_c.append(characters.remove(random_index));
            
           //..*..
        }
        
        //need to sort fixed_positions
        
        
        
         for(int index: fixed_positions.keySet())
         {
             char AA=fixed_positions.get(index);
             
             
             //System.out.println(index+".."+NON_SPECIFIC_AA);
             shuffled_c.insert(index, AA);
             
             
        
         }
        
                
        shuffled=shuffled_c.toString();
        System.out.println(sequence+"^^^^^^^^^^^^^^^^^^^^^^^^^^"+shuffled);
    }
    
    void calculate_pepmass()//check
    {
        //System.out.println(sequence);
        pepmass= deBruijn.WATER_MASS;
        
        //int length=sequence.length();
        for(int i=0; i<length; i++)
        {
            String s=sequence.substring(i, i+1);
            AminoAcid AA=(AminoAcid) deBruijn.AminoAcids.get(s);
            pepmass+=AA.residue_mass;
            
        }
        
    }
    
   

}
