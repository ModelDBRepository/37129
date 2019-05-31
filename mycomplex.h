typedef struct  
      {
        double realpart;
        double imgpart;
       }  mycomplex;

   mycomplex complexplus(mycomplex complex01, mycomplex complex02)
       {
         mycomplex complex_result;
         complex_result.realpart=complex01.realpart + complex02.realpart;
         complex_result.imgpart=complex01.imgpart + complex02.imgpart;
         return (complex_result);
       };

    mycomplex complexmulti(mycomplex complex01, mycomplex complex02)
       {
         mycomplex complex_result;
         complex_result.realpart=complex01.realpart * complex02.realpart
                                 -complex01.imgpart * complex02.imgpart;
         complex_result.imgpart=complex01.realpart * complex02.imgpart
                                +complex01.imgpart * complex02.realpart;
         return (complex_result);
       };

   mycomplex myconj(mycomplex complex01)
      {
       mycomplex complex_result;
       complex_result.realpart = complex01.realpart;
       complex_result.imgpart = -complex01.imgpart;
       return (complex_result);
      };
