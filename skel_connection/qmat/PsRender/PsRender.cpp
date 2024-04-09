//********************************************
// PsRender.h
// class CPsRenderer
//********************************************
// pierre.alliez@libertysurf.fr
// Adapted from :
//   Frederic Delhoume
//   Mark J. Kilgards 
// http://reality.sgi.com/opengl/tips/Feedback.html
//********************************************
// Created : 18/11/99
// Modified : 18/11/99
//********************************************

#define NOMINMAX

//#include "stdafx.h"
#include <math.h>
#include "PsRender.h"


//********************************************
// outputEPS
//********************************************
void CPsRenderer::Run(const char *pFilename,
											GLfloat *pFeedbackBuffer,
											int NbValues,
											BOOL sort)
{
//	ASSERT(pFilename);
  FILE *pFile = fopen(pFilename,"wt");
  if(!pFile)
	{
//		AfxMessageBox("Unable to open file for writing");
		return;
	}

  spewWireFrameEPS(pFile,sort,NbValues,pFeedbackBuffer);
  fclose(pFile);
}


// Code below is nothing more than Mark J. Kilgard's code
// adapted for C++/MFC

/* Copyright (c) Mark J. Kilgard, 1997. */

/* This program is freely distributable without licensing fees 
   and is provided without guarantee or warrantee expressed or 
   implied. This program is -not- in the public domain. */

/* Example showing how to use OpenGL's feedback mode to capture
   transformed vertices and output them as Encapsulated PostScript.
   Handles limited hidden surface removal by sorting and does
   smooth shading (albeit limited due to PostScript). */


//*********************************
// spewPrimitiveEPS
//*********************************
GLfloat *CPsRenderer::spewPrimitiveEPS(FILE * file, 
																			 GLfloat * loc)
{
  int token;
  int nvertices, i;
  GLfloat red, green, blue;
  int smooth;
  GLfloat dx, dy, dr, dg, db, absR, absG, absB, colormax;
  int steps;
  Feedback3Dcolor *vertex;
  GLfloat xstep, ystep, rstep, gstep, bstep;
  GLfloat xnext, ynext, rnext, gnext, bnext, distance;

  token = (int)*loc;
  loc++;
  switch (token) 
	{
  case GL_LINE_RESET_TOKEN:
  case GL_LINE_TOKEN:
    vertex = (Feedback3Dcolor *) loc;

    dr = vertex[1].red - vertex[0].red;
    dg = vertex[1].green - vertex[0].green;
    db = vertex[1].blue - vertex[0].blue;

    if(dr != 0 || dg != 0 || db != 0) 
		{
      /* Smooth shaded line. */
      dx = vertex[1].x - vertex[0].x;
      dy = vertex[1].y - vertex[0].y;
      distance = (float)sqrt(dx * dx + dy * dy);

      absR = (float)fabs(dr);
      absG = (float)fabs(dg);
      absB = (float)fabs(db);

      colormax = std::max(absR, std::max(absG, absB));
      steps = (int)std::max(1.0, colormax * distance * EPS_SMOOTH_LINE_FACTOR);

      xstep = dx / steps;
      ystep = dy / steps;

      rstep = dr / steps;
      gstep = dg / steps;
      bstep = db / steps;

      xnext = vertex[0].x;
      ynext = vertex[0].y;
      rnext = vertex[0].red;
      gnext = vertex[0].green;
      bnext = vertex[0].blue;

      /* Back up half a step; we want the end points to be
         exactly the their endpoint colors. */
      xnext -= (float)xstep / 2.0f;
      ynext -= (float)ystep / 2.0f;
      rnext -= (float)rstep / 2.0f;
      gnext -= (float)gstep / 2.0f;
      bnext -= (float)bstep / 2.0f;

			fprintf(file, "%g %g %g setrgbcolor\n",
				vertex[0].red, vertex[0].green, vertex[0].blue);
			fprintf(file, "%g %g moveto\n", vertex[0].x, vertex[0].y);

			for (i = 0; i < steps; i++) 
			{
				xnext += xstep;
				ynext += ystep;
				rnext += rstep;
				gnext += gstep;
				bnext += bstep;
				fprintf(file, "%g %g lineto stroke\n", xnext, ynext);
				fprintf(file, "%g %g %g setrgbcolor\n", rnext, gnext, bnext);
				fprintf(file, "%g %g moveto\n", xnext, ynext);
	    }
	    fprintf(file, "%g %g lineto stroke\n", vertex[1].x, vertex[1].y);
    } 
		else // single color segment -> S
		{
			fprintf(file, "%g %g %g %g %g %g %g S\n",vertex[1].x,vertex[1].y,
				vertex[0].x,vertex[0].y,vertex[0].red, vertex[0].green, vertex[0].blue);
		}

    loc += 14;          /* Each vertex element in the feedback
                           buffer is 7 GLfloats. */

    break;
  case GL_POLYGON_TOKEN:
    nvertices = (int)*loc;
    loc++;

    vertex = (Feedback3Dcolor *) loc;

    if(nvertices > 0) 
		{
      red = vertex[0].red;
      green = vertex[0].green;
      blue = vertex[0].blue;
      smooth = 0;
      for (i = 1; i < nvertices; i++) 
        if (red   != vertex[i].red || 
					  green != vertex[i].green || 
						blue  != vertex[i].blue) 
				{
          smooth = 1;
          break;
        }
      if(smooth) 
			{
        /* Smooth shaded polygon; varying colors at vertices. */
        int triOffset;

        /* Break polygon into "nvertices-2" triangle fans. */
        for (i = 0; i < nvertices - 2; i++) {
          triOffset = i * 7;

          fprintf(file, "[%g %g %g %g %g %g]",
            vertex[0].x, vertex[i + 1].x, vertex[i + 2].x,
            vertex[0].y, vertex[i + 1].y, vertex[i + 2].y);
          fprintf(file, " [%g %g %g] [%g %g %g] [%g %g %g] gouraudtriangle\n",
            vertex[0].red, vertex[0].green, vertex[0].blue,
            vertex[i + 1].red, vertex[i + 1].green, vertex[i + 1].blue,
            vertex[i + 2].red, vertex[i + 2].green, vertex[i + 2].blue);
					/*
					// x3 y3 r3 g3 b3 x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 ST
          fprintf(file,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ST\n",
						vertex[i+2].x,vertex[i+2].y,vertex[i+2].red,vertex[i+2].green,vertex[i+2].blue,
						vertex[i+1].x,vertex[i+1].y,vertex[i+1].red,vertex[i+1].green,vertex[i+1].blue,
						vertex[0].x,vertex[0].y,vertex[0].red,vertex[0].green,vertex[0].blue);*/
        }
      } 
			else 
			{
				if(nvertices > 3) // polygon
				{
					// Draw a filled polygon
					fprintf(file, "newpath\n");
					fprintf(file, "%g %g %g setrgbcolor\n", red, green, blue);
					fprintf(file, "%g %g moveto\n", vertex[0].x, vertex[0].y);
					for (i = 1; i < nvertices; i++)
						fprintf(file, "%g %g lineto\n", vertex[i].x, vertex[i].y);
					fprintf(file, "closepath fill\n\n");
				}
				else // triangle -> T
				{
					fprintf(file,"%g %g %g %g %g %g %g %g %g T\n",
						vertex[2].x,vertex[2].y,
						vertex[1].x,vertex[1].y,
						vertex[0].x,vertex[0].y,
						red,green,blue);
				}
      }
    }
    loc += nvertices * 7;  /* Each vertex element in the
                              feedback buffer is 7 GLfloats. */
    break;
  case GL_POINT_TOKEN:
    vertex = (Feedback3Dcolor *) loc;
    fprintf(file, "%g %g %g 0 360 %g %g %g P\n",vertex[0].x, vertex[0].y,
			m_PointSize / 2.0,vertex[0].red,vertex[0].green,vertex[0].blue);
    loc += 7;           /* Each vertex element in the feedback
                           buffer is 7 GLfloats. */
    break;
  default:
    /* XXX Left as an excersie to the reader. */
    printf("Incomplete implementation.  Unexpected token (%d).\n", token);
    exit(1);
  }
  return loc;
}

//*********************************
// spewUnsortedFeedback 
//*********************************
void CPsRenderer::spewUnsortedFeedback(FILE * file, 
																			 GLint size, 
																			 GLfloat * buffer)
{
  GLfloat *loc, *end;

  loc = buffer;
  end = buffer + size;
  while (loc < end) {
    loc = spewPrimitiveEPS(file, loc);
  }
}


//*********************************
// compare 
//*********************************
int CPsRenderer::compare(const void *a, 
												 const void *b)
{
  DepthIndex *p1 = (DepthIndex *) a;
  DepthIndex *p2 = (DepthIndex *) b;
  GLfloat diff = p2->depth - p1->depth;

  if (diff > 0.0) {
    return 1;
  } else if (diff < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

//*********************************
// spewSortedFeedback 
//*********************************
void CPsRenderer::spewSortedFeedback(FILE * file, 
																		 GLint size, 
																		 GLfloat * buffer)
{
  int token;
  GLfloat *loc, *end;
  Feedback3Dcolor *vertex;
  GLfloat depthSum;
  int nprimitives, item;
  DepthIndex *prims;
  int nvertices, i;

  end = buffer + size;

  /* Count how many primitives there are. */
  nprimitives = 0;
  loc = buffer;
  while (loc < end) {
    token = (int)*loc;
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      loc += 14;
      nprimitives++;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = (int)*loc;
      loc++;
      loc += (7 * nvertices);
      nprimitives++;
      break;
    case GL_POINT_TOKEN:
      loc += 7;
      nprimitives++;
      break;
    default:
      /* XXX Left as an excersie to the reader. */
      printf("Incomplete implementation.  Unexpected token (%d).\n",
        token);
      exit(1);
    }
  }

  /* Allocate an array of pointers that will point back at
     primitives in the feedback buffer.  There will be one
     entry per primitive.  This array is also where we keep the
     primitive's average depth.  There is one entry per
     primitive  in the feedback buffer. */
  prims = (DepthIndex *) malloc(sizeof(DepthIndex) * nprimitives);

  item = 0;
  loc = buffer;
  while (loc < end) {
    prims[item].ptr = loc;  /* Save this primitive's location. */
    token = (int)*loc;
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z + vertex[1].z;
      prims[item].depth = (float)depthSum / 2.0f;
      loc += 14;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = (int)*loc;
      loc++;
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z;
      for (i = 1; i < nvertices; i++) {
        depthSum += vertex[i].z;
      }
      prims[item].depth = depthSum / nvertices;
      loc += (7 * nvertices);
      break;
    case GL_POINT_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      prims[item].depth = vertex[0].z;
      loc += 7;
      break;
    default:
      /* XXX Left as an excersie to the reader. */
//      ASSERT(1);
		;
    }
    item++;
  }
//  ASSERT(item == nprimitives);

  /* Sort the primitives back to front. */
  qsort(prims, nprimitives, sizeof(DepthIndex), compare);

  /* XXX Understand that sorting by a primitives average depth
     doesn't allow us to disambiguate some cases like self
     intersecting polygons.  Handling these cases would require
     breaking up the primitives.  That's too involved for this
     example.  Sorting by depth is good enough for lots of
     applications. */

  /* Emit the Encapsulated PostScript for the primitives in
     back to front order. */
  for (item = 0; item < nprimitives; item++) {
    (void) spewPrimitiveEPS(file, prims[item].ptr);
  }

  free(prims);
}


//*********************************
// spewWireFrameEPS 
//*********************************
void CPsRenderer::spewWireFrameEPS(FILE * file, 
																	 BOOL  doSort, 
																	 GLint size, 
																	 GLfloat * buffer)
{
  GLfloat clearColor[4], viewport[4];
  GLfloat lineWidth;

  // Read back a bunch of OpenGL state to help make the EPS
  // consistent with the OpenGL clear color, line width, point
  // size, and viewport.
  glGetFloatv(GL_VIEWPORT, viewport);
  glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColor);
  glGetFloatv(GL_LINE_WIDTH, &lineWidth);
  glGetFloatv(GL_POINT_SIZE, &m_PointSize);

  /* Emit EPS header. */
  fputs("%!PS-Adobe-2.0 EPSF-2.0\n", file);
  /* Notice %% for a single % in the fprintf calls. */
  fprintf(file, "%%%%BoundingBox: %g %g %g %g\n",
    viewport[0], viewport[1], viewport[2], viewport[3]);
  fputs("%%EndComments\n", file);
  fputs("\n", file);
  fputs("gsave\n", file);
  fputs("\n", file);

  fprintf(file, "/threshold %g def\n", EPS_GOURAUD_THRESHOLD);

  fprintf(file,"\n%% RGB color command - r g b C\n");
  fprintf(file,"/C { setrgbcolor } bind def\n");

	// This is more compact
  fprintf(file,"\n%% Point - x_center y_center PointSize/2 0 360 r g b P\n");
  fprintf(file,"/P { C arc fill } bind def\n");

  fprintf(file,"\n%% Segment - x2 y2 x1 y1 r g b S\n");
  fprintf(file,"/S { C moveto lineto stroke } bind def\n");

  fprintf(file,"\n%% Flat-shaded triangle - x3 y3 x2 y2 x1 y1 r g b T\n");
  fprintf(file,"/T { C newpath moveto lineto lineto closepath fill } bind def\n");

  // Output Frederic Delhoume's "gouraudtriangle" PostScript fragment.
  fputs("% The gouraudtriangle PostScript fragement below is free\n",file);
  fputs("% written by Frederic Delhoume (delhoume@ilog.fr)\n",file);

  for(int i = 0;gouraudtriangleEPS[i]; i++) 
    fprintf(file, "%s\n", gouraudtriangleEPS[i]);

  fprintf(file, "\n%g setlinewidth\n", EPS_LINE_WIDTH);

  // Clear the background like OpenGL had it
  fprintf(file, "%g %g %g setrgbcolor\n",
    clearColor[0], clearColor[1], clearColor[2]);
  fprintf(file, "%g %g %g %g rectfill\n\n",
    viewport[0], viewport[1], viewport[2], viewport[3]);

	// Main process
  if(doSort)
    spewSortedFeedback(file, size, buffer);
  else
    spewUnsortedFeedback(file, size, buffer);

  // Emit EPS trailer
  fputs("\ngrestore\n", file);
  fputs("showpage\n",file);

}

// ** EOF **
