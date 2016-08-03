



typedef struct                      /* Polygon vertex structure          */
{
  double              x;            /* Vertex x component                */
  double              y;            /* vertex y component                */
} Vertex;

typedef struct                      /* Vertex list structure             */
{
  int                 num_vertices; /* Number of vertices in list        */
  Vertex         *vertex;       /* Vertex array pointer              */
} VertexList;

typedef struct                      /* Polygon set structure             */
{
  VertexList    *contour;      /* Contour array pointer             */
} Polygon;





extern int resample_img(void);

extern void super_poly_clip( Polygon *ref_poly, Polygon *tmp_clip_poly, Polygon *clip_poly, long mm, long nn ) ;



