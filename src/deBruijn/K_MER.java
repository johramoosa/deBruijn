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
public class K_MER {
    
    Character decoy;/* replacement AA */
    int outgoing;
    
    K_MER()
    {
        outgoing=0;
    }
    
    K_MER(Character aa)
    {
        decoy=aa;
        outgoing=1;
    }
    
    K_MER(Character aa, int count)
    {
        decoy=aa;
        outgoing=count;
    }
    
    int add()
    {
        outgoing++;
        return outgoing;
    }
    
}
