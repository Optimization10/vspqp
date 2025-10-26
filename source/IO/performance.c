#include "shared_api.h"


void VSP_Print_Performance()
{
    #if MEASURE_PERFORMANCE
        FILE *file = fopen("output/performance.csv", "w") ;
        fprintf(file, ",Time,Count") ;
        fprintf(file,"\nMC,%Lf,%llu", G_Performance_Struct.MC_Seconds, G_Performance_Struct.MC_Count) ;
        fprintf(file,"\nFM,%Lf,%llu", G_Performance_Struct.FM_Seconds, G_Performance_Struct.FM_Count) ;
        fprintf(file,"\nFlips,%Lf,%llu",G_Performance_Struct.Flips_Seconds, G_Performance_Struct.Flips_Count) ;
        fprintf(file,"\nAll Coarsening,%Lf,%llu", G_Performance_Struct.All_Coarsening_Seconds, G_Performance_Struct.All_Coarsening_Count) ;
        fprintf(file,"\nFM_Flips,%Lf, %llu", G_Performance_Struct.FM_Flips_Seconds, G_Performance_Struct.FM_Flips_Count) ;
        fprintf(file,"\nRelaxed_To_VSP,%Lf,%llu", G_Performance_Struct.Relaxed_To_VSP_Seconds, G_Performance_Struct.Relaxed_To_VSP_Count) ;

        fclose(file) ;
    #endif
    return ;
}