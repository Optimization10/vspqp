#ifndef SHARED_MACROS_H
#define SHARED_MACROS_H






#define NO_FLIP -1
#define XUP 0
#define YUP 1
#define XDN 2
#define YDN 3



#define NONE -1
#define X0 0
#define Y0 1
#define X1 2
#define Y1 3

#define NO_VAR -1
#define X_VAR   0
#define Y_VAR   1

// Given a variable 0 or 1, and its value in an entry 0 or 1, get the encoding
// associated that entry
// var = 0 = x, zi = 0 --> 0 X0
// var = 0 = x, zi = 1 --> 2 X1
// var = 1 = y, zi = 0 --> 1 Y0
// var = 1 = y, zi = 1 --> 3 Y1

#define UNLOCKED(locked, i) (  ((locked) == NULL) || ( !(locked)[(i)] )  )

// Swap the values of A and B
#define SWAP(a, b, type) do {                     \
    type VSPTempForSwappingONLYDoNOTUse = (a);    \
    (a) = (b) ;                                   \
    (b) = VSPTempForSwappingONLYDoNOTUse ;        \
} while (0)


#define IMBALANCE(wgtz, wgtw) \
({ \
    typeof(wgtz) _diff = (wgtz) - (wgtw); \
    _diff < 0 ? -_diff : _diff; \
})

 

// A, B are logical conditions
#define IFF(A, B) ( ((A) + (B) == 0) || ((A) + (B) == 2) )

// Here, x is one of X0, Y0, X1, Y1
// X0, X1 (that is 0 and 2) map to 0 (X_VAR), and Y0 and Y1 (that is 1 and 3) map to 1 (Y_VAR) 
#define GET_VAR(x) ( (x) % 2)



// Plug in 0, 1, 2, or 3 according to our encoding
// X0 ---> 0
// Y0 ---> 0
// X1 ---> 1
// Y1 ---> 1
#define GET_ZI(x) ( x >= 2 )

// X0 = 0 --> 1 = Y_VAR
// Y0 = 1 --> 0 = X_VAR
// X1 = 2 --> 1 = Y_VAR
// Y1 = 3 --> 0 = X_VAR
#define GET_OTHER_VAR(flip)   ( ( (flip) + 1) % 2 )
#define IS_UPFLIP(flip) ( (flip) <= 1 )
#define IS_DNFLIP(flip) ( (flip) >= 2 )

// 0 --> 1
// 1 --> 0
#define OTHER_VAR(v) ( 1 - (v) )

// Get Aw[i], where w is the variable NOT associated with "flip"
#define GET_AWI(flip, i, Avars) ( (Avars)[GET_OTHER_VAR( (flip) )][ (i) ] )

#define DBG_INFO __FILE__, __LINE__, count
#define ERR(msg) Error( msg, __FILE__, __LINE__, count)
#define CAST(arg, type) ( ((type) (arg)) )

// Plug in 0, 1, 2, or 3 according to our encoding
// X0 --->  1
// Y0 --->  1
// X1 ---> -1
// Y1 ---> -1
#define GET_DIR(x) ( ( (x) < 2 ) + (x >= 2)*(-1) )

// X0 ---> X1
// Y0 ---> Y1
// X1 ---> X0
// Y1 ---> Y0
#define FLIP_COMP(x) ( ( (x) + 2) % 4)

#define FEASIBLE_UPFLIP(wgt, ci, u) ( (wgt) + (ci) <= (u) )

// Lower bound of 1 is disabled
#define FEASIBLE_DNFLIP(wgt, ci) ( 1 )

#define IS_FEASIBLE_FLIP(zi, wgt, ci, u) ( zi ? FEASIBLE_DNFLIP(wgt, ci) : FEASIBLE_UPFLIP(wgt, ci, u) )

// zi = 0 -->  1
// zi = 1 --> -1
#define DIR_FLIP(zi) ( 1 - 2 * (zi) )


#define UNUSED __attribute__((unused))

#define GAIN_FLIP(dir, ci, gamma, Awi) ( (dir) * ( (ci) - (gamma) * (Awi) ) )

#ifdef DEBUG
    #define DBGCALL(f) (f)
#else
    #define DBGCALL(f)
#endif

/* Return true if low <= num < high, and false otherwise */
#define IN_RANGE(num, low, high) ( ( (low) <= (num) ) && ( (num) < (high) ) ) 

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/* Returns true iff a/b < c/d */
#define CMP_FRAC(a, b, c, d) ((a) * (d) < (c) * (b))
#define CMPEQ_FRAC(a, b, c, d) ((a) * (d) <= (c) * (b))

#ifdef DEBUG
    #define VSP_ASSERT( cond ) \
        do { \
        if ( !(cond) ) { \
           Error("Failed assert.", __FILE__, __LINE__, count) ; \
        } \
    } while (0)
#else
    #define VSP_ASSERT( cond ) 
#endif

#if defined(DEBUG) || defined(DEBUG_COARSENING)
    #define WARN_IF_NOT_DEBUG_MODE
#else
    #define WARN_IF_NOT_DEBUG_MODE printf("\n\nWARNING: Calling debugging function when not in debug mode. File: %s  Line: %lld\n\n", __FILE__, (long long) __LINE__)
#endif



#ifdef DEBUG
    #define SET_FUNCTION_COUNTER static long long unsigned count = 0 ; count++
#else
    #define SET_FUNCTION_COUNTER
#endif


#endif