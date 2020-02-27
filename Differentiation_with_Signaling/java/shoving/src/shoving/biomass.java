package shoving;

import static java.lang.Math.*;
//import java.util.Vector;

public class biomass {

  public Results pushing3D(int nbac, double[] bac_x, double[] bac_y, double[] bac_z, double[] bac_r) {

        int push = 1;       // boolean test for more shoving (1: shoving needed; 0: shoving ready)
        double[] shx = new double[nbac];
        double[] shy = new double[nbac]; 
        double[] shz = new double[nbac]; 

        double dx, dy,dz, d, overlap;

        while (push == 1) { // One shoving step
            push = 0;       // start assuming there is no need for pushing

            for (int i = 0; i < nbac; i++) {   // browse n cells
                for (int j = i + 1; j < nbac; j++) {  // search overlapping cells among all other cells ...
                    dx = bac_x[i] - bac_x[j];
                    dy = bac_y[i] - bac_y[j]; 
                    dz = bac_z[i] - bac_z[j]; 
                    // calculate euclidian distance d between cell i (reference cell) and the new cell j
                    d = sqrt(dx*dx + dy*dy+ dz*dz );
                    // calculate overlap
                    overlap = bac_r[i]+bac_r[j]-d;
                                
                    if (overlap>1e-18){
                        shx[i] = shx[i] + dx*overlap/d;
                        shy[i] = shy[i] + dy*overlap/d; 
                        shz[i] = shz[i] + dz*overlap/d; 
                    
                        shx[j] = shx[j] - dx*overlap/d;
                        shy[j] = shy[j] - dy*overlap/d; 
                        shz[j] = shz[j] - dz*overlap/d; 
                        push = 1;
                    }
                }
            }

            for (int i = 0; i < nbac; i++) {
                // move cell i with the resultant components in x and y directions
                bac_x[i] = bac_x[i] + shx[i];
                bac_y[i] = bac_y[i] + shy[i]; 
                bac_z[i] = bac_z[i] + shz[i]; 

                

                // reset the resultant vectors of movement to zero
                shx[i] = 0;
                shy[i] = 0; 
                shz[i] = 0; 
            }
        }

        Results Results = new Results();
        Results.bac_x = bac_x;
        Results.bac_y = bac_y; 
        Results.bac_z = bac_z; 
        return Results;
    }     

}
