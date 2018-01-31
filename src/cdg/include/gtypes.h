//////////////////////////////////////////////////////////////
// MODULE   : gtypes 
// DESC     : basic data types
// DATE     : 1/1/2018
// AUTHOR   : DR
// COPYRIGHT: CIRA/CSU 2018-2020
////////////////////////////////////////////////////////////////
#include <cstddef>


// Globally-defined defs:
// _G_AUTO_CREATE_DEV: Auto-create class data and class on device if ACC defined
// _G_BOUNDS_CHK     : Do data bounds check

// Basic data types:
#define GCHAR      char
#define GSHORT     short
#define GUSHORT    unsigned short
#define GINT       int
#define GBOOL      bool
#define GUINT      unsigned int
#define GLONG      long
#define GLONGLONG  long long
#define GSIZET     size_t
#define GFLOAT     float
#define GDOUBLE    double
#define GQUAD      long double


#if !defined(GBOOL)
#define GBOOL bool
#endif
#if !defined(TRUE)
#define TRUE  true
#endif
#if !defined(FALSE)
#define FALSE false
#endif

// Miscellaneous defs:
#define NULLPTR    NULL
#define PI         3.14159265358979323846264338328

// Datatypes for communication:
#if !defined(GC_DATATYPE_DEFTYPE_DEF)
#define  GC_DATATYPE_DEFTYPE_DEF
enum GC_DATATYPE_DEFTYPE {G_GFLOAT=0,G_GDOUBLE ,G_GQUAD  ,G_GINT  ,G_GSHORT ,G_GUSHORT ,
                          G_GLONG   ,G_BYTE    ,G_GUCHAR ,G_GWORD ,G_GDWORD ,G_GFPOS   ,
                          G_GKEY    ,G_GNODEID ,G_ELEMTYPE,G_BDYTYPE};
const GINT  G_TYPESZ[] = {sizeof  (GFLOAT),sizeof  (GDOUBLE),sizeof (GQUAD),sizeof (GINT ),sizeof(GSHORT ),
                          sizeof(GUSHORT ),sizeof   (GLONG ),sizeof (GBYTE),sizeof(GUCHAR),
                          sizeof   (GWORD),sizeof   (GDWORD),sizeof (GFPOS),sizeof  (GKEY),
                          sizeof (GNODEID),sizeof    (GINT ),sizeof (GINT )};
#endif

