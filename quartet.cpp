//
//  quartet.cpp
//  iqtree
//
//  Created by Minh Bui on 24/07/15.
//
//

#include <stdio.h>

#include "phylotree.h"
#include "phylosupertree.h"
#include "phylosupertree.h"
// #include "lmap.c"

#if 0  /*** moved to phylotree.h ***/
/* Index definition for counter array needed in likelihood mapping analysis (HAS) */
#define LM_REG1 0   /* top corner */
#define LM_REG2 1   /* bottom-right corner */
#define LM_REG3 2   /* bottom-left corner */
#define LM_REG4 3   /* right rectangle */
#define LM_REG5 4   /* bottom rectangle */
#define LM_REG6 5   /* left rectangle */
#define LM_REG7 6   /* center */
#define LM_AR1  7   /* top third */
#define LM_AR2  8   /* bottom-right third */
#define LM_AR3  9   /* bottom-left third */
#define LM_MAX  10
#endif

//*** likelihood mapping stuff (imported from TREE-PUZZLE's lmap.c) (HAS)

// #include <time.h>

/**********************************************************/
/*  Likelihood mapping routines (TODO: move to lmap.c/h)  */
/**********************************************************/

/*
                   (a,b)-(c,d)   => numclust == 4
                   (a,b)-(c,c)   => numclust == 3
                   (a,a)-(b,b)   => numclust == 2

                     1l/\1r
                      /  \
                     /  1 \
                 6u / \  / \ 4u
                   /   \/   \
               6d /    /\    \ 4d
                 / 6  /  \  4 \
                /\   /  7 \   /\
             3u/  \ /______\ /  \2u
              / 3  |    5   |  2 \
             /_____|________|_____\
	        3d   5l  5r   2d
    (a,d)-(b,c)                  (a,c)-(b,d)   => numclust == 4
    (a,c)-(b,c)                  (a,c)-(b,c)   => numclust == 3
    (a,b)-(a,b)                  (a,b)-(a,b)   => numclust == 2

*/

/***********************************
*  Likelihood mapping to SVG file  *
***********************************/

/* first lines of SVG likelihood mapping file */
void initsvg(FILE *ofp)
{
	/* SVG preamble */
	fprintf(ofp,"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
	fprintf(ofp,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
	fprintf(ofp,"<svg\n");
	fprintf(ofp,"   xmlns:svg=\"http://www.w3.org/2000/svg\"\n");
	fprintf(ofp,"   xmlns=\"http://www.w3.org/2000/svg\"\n");
	fprintf(ofp,"   xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n");
	fprintf(ofp,"   version=\"1.1\"\n");
	fprintf(ofp,"   baseProfile=\"full\"\n");
	fprintf(ofp,"   id=\"body\"\n");
	fprintf(ofp,"   width=\"800px\"\n");
	fprintf(ofp,"   height=\"800px\"\n");
	fprintf(ofp,"   viewBox=\"0 0 1000 1000\"\n");
	fprintf(ofp,"   preserveAspectRatio=\"none\">\n");
	fprintf(ofp,"  <defs>\n");
	fprintf(ofp,"    <style type=\"text/css\"><![CDATA[\n");
	fprintf(ofp,"      circle{ stroke: none; }\n");
	fprintf(ofp,"      polygon{ stroke: black; stroke-width: 2px; fill: none; }\n");
	fprintf(ofp,"      line{ stroke: black; stroke-width: 2px; }\n");
	fprintf(ofp,"      text{ font-size:50px; }\n");
	fprintf(ofp,"    ]]></style>\n");
	fprintf(ofp,"  </defs>\n");
	fprintf(ofp,"  <title\n");
	fprintf(ofp,"     id=\"title1\">SVG drawing</title>\n");
	/* end SVG preamble */

	/* triangle 1 (top) */
	fprintf(ofp,"<g transform=\"scale(0.45)\"><g transform=\"translate(600,1050)\">\n");
	fprintf(ofp,"  <g id=\"fig1\">\n");
	fprintf(ofp,"	<polygon points=\"0.0,-0.0 1000.0,-0.0 500,-866.0254038\" />\n");


#if LMAP_CLUSTER
	if (numclust == 4) { /* four cluster analysis */
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"500.0\"\n");
		fprintf(ofp,"	   y=\"-896.0254038\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_top_1\">(a,b)-(c,d)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"-30.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_left_1\">(a,d)-(b,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"1030.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_right_1\">(a,c)-(b,d)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
	}
	if (numclust == 3) { /* three cluster analysis */
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"500.0\"\n");
		fprintf(ofp,"	   y=\"-896.0254038\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_top_1\">(a,b)-(c,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"-30.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_left_1\">(a,c)-(b,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"1030.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_right_1\">(a,c)-(b,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
	}
	if (numclust == 2) { /* two cluster analysis */
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"500.0\"\n");
		fprintf(ofp,"	   y=\"-896.0254038\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_top_1\">(a,a)-(b,b)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"-30.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_left_1\">(a,b)-(a,b)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"1030.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		fprintf(ofp,"	   id=\"label_right_1\">(a,b)-(a,b)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
	}
#endif /* LMAP_CLUSTER */

} /* initsvg */




void plotlmpointsvg(FILE *ofp, double w1, double w2)
{
	/* plot dots into triangle 1 (top) */
	fprintf(ofp,"	<circle cx=\"%.10f\" cy=\"%.10f\" r=\"2\" />\n", (0.5*w1 + w2)*1000, -(w1*866.0254038));
} /* plotlmpointsvg */



// void finishsvg(FILE *ofp, unsigned long **countarr)
void finishsvg(FILE *ofp, vector<SeqQuartetInfo> lmap_seq_quartet_info, int leafNum, unsigned long Numquartets)
{
	fprintf(ofp,"  </g>\n");
	/* end triangle 1 (top) */

	/* triangle 2 (bottom left) */
	fprintf(ofp,"  <g id=\"fig2\" transform=\"translate(-550.0,1000)\">\n");
	fprintf(ofp,"	<polygon points=\"0.0,-0.0 1000.0,-0.0 500.0,-866.0254038\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line2-1\"\n");
	fprintf(ofp,"	   y2=\"-0.0\"\n");
	fprintf(ofp,"	   x2=\"500\"\n");
	fprintf(ofp,"	   y1=\"-288.6751346\"\n");
	fprintf(ofp,"	   x1=\"500\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line2-2\"\n");
	fprintf(ofp,"	   y2=\"-433.0127019\"\n");
	fprintf(ofp,"	   x2=\"250\"\n");
	fprintf(ofp,"	   y1=\"-288.6751346\"\n");
	fprintf(ofp,"	   x1=\"500\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line2-3\"\n");
	fprintf(ofp,"	   y2=\"-433.0127019\"\n");
	fprintf(ofp,"	   x2=\"750\"\n");
	fprintf(ofp,"	   y1=\"-288.6751346\"\n");
	fprintf(ofp,"	   x1=\"500\" />\n");
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"440\"\n");
	fprintf(ofp,"	   y=\"-500\"\n");
	fprintf(ofp,"	   id=\"up_2\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_AR1]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"up_2\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_AR1]*100.0/Numquartets);
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"250\"\n");
	fprintf(ofp,"	   y=\"-150\"\n");
	fprintf(ofp,"	   id=\"down_left_2\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_AR3]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"down_left_2\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_AR3]*100.0/Numquartets);
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"630\"\n");
	fprintf(ofp,"	   y=\"-150\"\n");
	fprintf(ofp,"	   id=\"down_right_2\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_AR2]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"down_right_2\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_AR2]*100.0/Numquartets);
	fprintf(ofp,"  </g>\n");
	/* end triangle 2 (bottom left) */

	/* triangle 3 (bottom right) */
	fprintf(ofp,"  <g id=\"fig3\" transform=\"translate(550,1000)\">\n");
	fprintf(ofp,"	<polygon points=\"0.0,-0.0 1000.0,-0.0 500.0,-866.0254038\" />\n");
	fprintf(ofp,"	<polygon id=\"triangle3b\" points=\"250,-144.3375673 750,-144.3375673 500,-577.3502692\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line3-1\"\n");
	fprintf(ofp,"	   x1=\"125\"\n");
	fprintf(ofp,"	   y1=\"-216.5063509\"\n");
	fprintf(ofp,"	   x2=\"250\"\n");
	fprintf(ofp,"	   y2=\"-144.3375673\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line3-2\"\n");
	fprintf(ofp,"	   x1=\"375\"\n");
	fprintf(ofp,"	   y1=\"-649.5190528\"\n");
	fprintf(ofp,"	   x2=\"500\"\n");
	fprintf(ofp,"	   y2=\"-577.3502692\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line3-3\"\n");
	fprintf(ofp,"	   x1=\"625\"\n");
	fprintf(ofp,"	   y1=\"-649.5190528\"\n");
	fprintf(ofp,"	   x2=\"500\"\n");
	fprintf(ofp,"	   y2=\"-577.3502692\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line3-4\"\n");
	fprintf(ofp,"	   x1=\"875\"\n");
	fprintf(ofp,"	   y1=\"-216.5063509\"\n");
	fprintf(ofp,"	   x2=\"750\"\n");
	fprintf(ofp,"	   y2=\"-144.3375673\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line3-5\"\n");
	fprintf(ofp,"	   x1=\"750\"\n");
	fprintf(ofp,"	   y1=\"-0.0\"\n");
	fprintf(ofp,"	   x2=\"750\"\n");
	fprintf(ofp,"	   y2=\"-144.3375673\" />\n");
	fprintf(ofp,"	<line\n");
	fprintf(ofp,"	   id=\"line3-6\"\n");
	fprintf(ofp,"	   x1=\"250\"\n");
	fprintf(ofp,"	   y1=\"-0.0\"\n");
	fprintf(ofp,"	   x2=\"250\"\n");
	fprintf(ofp,"	   y2=\"-144.3375673\" />\n");

	/* number of resolved quartets, top */
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"500\"\n");
	fprintf(ofp,"	   y=\"-660\"\n");
	fprintf(ofp,"	   text-anchor=\"middle\"\n");
	fprintf(ofp,"	   id=\"up_3\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_REG1]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"up_3\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_REG1]*100.0/Numquartets);

	/* number of resolved quartets, bottom left */
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   y=\"-50\"\n");
	fprintf(ofp,"	   x=\"70\"\n");
	fprintf(ofp,"	   id=\"down_left_3\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_REG3]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"down_left_3\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_REG3]*100.0/Numquartets);

	/* number of resolved quartets, bottom right */
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   y=\"-50\"\n");
	fprintf(ofp,"	   x=\"770\"\n");
	fprintf(ofp,"	   id=\"down_right_3\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_REG2]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"down_right_3\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_REG2]*100.0/Numquartets);

	/* number of partly resolved quartets, bottom */
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"500\"\n");
	fprintf(ofp,"	   y=\"-50\"\n");
	fprintf(ofp,"	   text-anchor=\"middle\"\n");
	fprintf(ofp,"	   id=\"down_side_3\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_REG5]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"down_side_3\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_REG5]*100.0/Numquartets);

	/* number of unresolved quartets, center */
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"500\"\n");
	fprintf(ofp,"	   y=\"-280\"\n");
	fprintf(ofp,"	   text-anchor=\"middle\"\n");
	fprintf(ofp,"	   id=\"center_3\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_REG7]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"center_3\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_REG7]*100.0/Numquartets);

	/* number of partly resolved quartets, top right */
	/* fprintf(ofp,"<circle cx=\"685.0\" cy=\"-390.8439\" r=\"20\" />\n"); */ /* ro */
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"685.0\"\n");
	fprintf(ofp,"	   y=\"-390.8439\"\n");
	fprintf(ofp,"	   text-anchor=\"middle\"\n");
	fprintf(ofp,"	   transform=\"rotate(60,665.0,-380.8439)\"\n");
	fprintf(ofp,"	   id=\"right_side_3\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_REG4]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"right_side_3\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_REG4]*100.0/Numquartets);

	/* number of partly resolved quartets, top left */
	/* fprintf(ofp,"<circle cx=\"315.0\" cy=\"-390.8439\" r=\"20\" />\n"); */  /* lo */
	fprintf(ofp,"	<text\n");
	fprintf(ofp,"	   x=\"315.0\"\n");
	fprintf(ofp,"	   y=\"-390.8439\"\n");
	fprintf(ofp,"	   text-anchor=\"middle\"\n");
	fprintf(ofp,"	   transform=\"rotate(-60,335.0,-380.8439)\"\n");
	fprintf(ofp,"	   id=\"left_side_3\">%.1f%%</text>\n", (double)lmap_seq_quartet_info[leafNum].countarr[LM_REG6]*100.0/Numquartets);
	// fprintf(ofp,"	   id=\"left_side_3\">%.1f%%</text>\n", (double)countarr[Maxspc][LM_REG6]*100.0/Numquartets);

	fprintf(ofp,"  </g>\n");
	/* end triangle 3 (bottom right) */

	fprintf(ofp,"</g></g>\n");
	fprintf(ofp,"</svg>\n");
} /* finishsvg */

/* end - Likelihood mapping to SVG file  */


/***********************************
*  Likelihood mapping to EPS file  *
***********************************/

/* first lines of EPSF likelihood mapping file */
void initeps(FILE *ofp)
{
	time_t Starttime;
	time(&Starttime);

	fprintf(ofp, "%%!PS-Adobe-3.0 EPSF-3.0\n");
	fprintf(ofp, "%%%%BoundingBox: 60 210 550 650\n");
	fprintf(ofp, "%%%%Pages: 1\n");
	fprintf(ofp, "%%%%Creator: IQ-TREE/TREE-PUZZLE\n");
#if 0
#	ifndef ALPHA
		fprintf(ofp, "%%%%Creator: %s (version %s)\n", PACKAGE, VERSION);
#	else
		fprintf(ofp, "%%%%Creator: %s (version %s%s)\n", PACKAGE, VERSION, ALPHA);
#	endif
#endif
	fprintf(ofp, "%%%%Title: Likelihood Mapping Analysis\n");
	fprintf(ofp, "%%%%CreationDate: %s", asctime(localtime(&Starttime)) );
	fprintf(ofp, "%%%%DocumentFonts: Helvetica\n");
	fprintf(ofp, "%%%%DocumentNeededFonts: Helvetica\n");
	fprintf(ofp, "%%%%EndComments\n");
	fprintf(ofp, "%% use inch as unit\n");
	fprintf(ofp, "/inch {72 mul} def\n");
	fprintf(ofp, "%% triangle side length (3 inch)\n");
	fprintf(ofp, "/tl {3 inch mul} def\n");
	fprintf(ofp, "%% plot one dot (x-y coordinates on stack)\n");
	fprintf(ofp, "/dot {\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "0.002 tl 0 360 arc  %% radius is 0.002 of the triangle length\n");
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "fill\n");
	fprintf(ofp, "} def\n");

	/* PS definition of a flush right print */
	fprintf(ofp, "\n%% flush right show\n");
	fprintf(ofp, "/centershow {\n");
	fprintf(ofp, "   dup stringwidth pop  %% get length of string\n");
	fprintf(ofp, "   neg 0 rmoveto        %% move width to left\n");
	fprintf(ofp, "   show\n");
	fprintf(ofp, "} def\n");
	fprintf(ofp, "\n%% centered show\n");

	/* PS definition of a centered print */
	fprintf(ofp, "/centershow {\n");
	fprintf(ofp, "   dup stringwidth pop %% get length of string\n");
	fprintf(ofp, "   -2 div              %% devide length by -2\n");
	fprintf(ofp, "   0 rmoveto           %% move half width to left\n");
	fprintf(ofp, "   show\n");
	fprintf(ofp, "} def\n");


	fprintf(ofp, "%% preamble\n");
	fprintf(ofp, "/Helvetica findfont\n");
	fprintf(ofp, "12 scalefont\n");
	fprintf(ofp, "setfont\n");
	fprintf(ofp, "%% 0/0 for triangle of triangles\n");
	fprintf(ofp, "0.9 inch 3 inch translate\n");
	fprintf(ofp, "%% first triangle (the one with dots)\n");
	fprintf(ofp, "0.6 tl 1.2 tl 0.8660254038 mul translate\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.0 tl 0.0 tl moveto\n");
	fprintf(ofp, " 1.0 tl 0.0 tl lineto\n");
	fprintf(ofp, " 0.5 tl 0.8660254038 tl lineto\n");
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");

#if LMAP_CLUSTER
	if (numclust == 4) { /* four cluster analysis */
		fprintf(ofp, "%% label corners\n");
		fprintf(ofp, "0.5 tl 0.9 tl moveto\n"); /* old: 0.375 0.9 */
		fprintf(ofp, "((a,b)-(c,d)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "-0.045 tl -0.08 tl moveto\n"); /* old: -0.16 -0.08 */
		fprintf(ofp, "((a,d)-(b,c)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "1.045 tl -0.08 tl moveto\n"); /* old: -0.92 -0.08 */
		fprintf(ofp, "((a,c)-(b,d)) centershow %% CHANGE HERE IF NECESSARY\n");
	}
	if (numclust == 3) { /* three cluster analysis */
		fprintf(ofp, "%% label corners\n");
		fprintf(ofp, "0.5 tl 0.9 tl moveto\n"); /* old: 0.375 0.9 */
		fprintf(ofp, "((a,b)-(c,c)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "-0.045 tl -0.08 tl moveto\n"); /* old: -0.16 -0.08 */
		fprintf(ofp, "((a,c)-(b,c)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "1.045 tl -0.08 tl moveto\n"); /* old: -0.92 -0.08 */
		fprintf(ofp, "((a,c)-(b,c)) centershow %% CHANGE HERE IF NECESSARY\n");
	}
	if (numclust == 2) { /* two cluster analysis */
		fprintf(ofp, "%% label corners\n");
		fprintf(ofp, "0.5 tl 0.9 tl moveto\n"); /* old: 0.375 0.9 */
		fprintf(ofp, "((a,a)-(b,b)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "-0.045 tl -0.08 tl moveto\n"); /* old: -0.16 -0.08 */
		fprintf(ofp, "((a,b)-(a,b)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "1.045 tl -0.08 tl moveto\n"); /* old: -0.92 -0.08 */
		fprintf(ofp, "((a,b)-(a,b)) centershow %% CHANGE HERE IF NECESSARY\n");
	}
#endif /* LMAP_CLUSTER */

} /* initeps */

/* plot one point of likelihood mapping analysis (EPS) */
void plotlmpointeps(FILE *epsofp, double w1, double w2)
{
	fprintf(epsofp,"%.10f tl %.10f tl dot\n", 0.5*w1 + w2, w1*0.8660254038);
} /* plotlmpointeps */



#if 0
/* plot one point of likelihood mapping analysis */
void plotlmpoint(FILE *epsofp, FILE *svgofp, double w1, double w2)
{
	if (lmapeps_optn) {
		fprintf(epsofp,"%.10f tl %.10f tl dot\n",
			0.5*w1 + w2, w1*0.8660254038);
	}
	if (lmapsvg_optn) {
		//fprintf(svgofp,"    <use x=\"%.10f\" y=\"%.10f\" xlink:href=\"#result\" />\n", (0.5*w1 + w2), -(w1*0.8660254038));
		fprintf(svgofp,"    <circle cx=\"%.10f\" cy=\"%.10f\" r=\"2\" />\n", 
			(0.5*w1 + w2)*1000, -(w1*866.0254038));
	}
} /* plotlmpoint */
#endif



#if 0
/* plot one point of likelihood mapping analysis */
void plotlmpointcolor(FILE *epsofp, FILE *svgofp, double w1, double w2, int red, int green, int blue)
{
	if (lmapeps_optn) {
		fprintf(epsofp,"currentrgbcolor %d %d %d setrgbcolor\n", red, green, blue);
		fprintf(epsofp,"%.10f tl %.10f tl dot\n",
			0.5*w1 + w2, w1*0.8660254038);
		fprintf(epsofp,"setrgbcolor\n");
	}
	if (lmapsvg_optn) {
/*
	stijn imbrechts:
		Adding colour to elements is pretty easy, if you are familiar with
		CSS, it works almost exactly the same.
		The dots are represented by a <circle> element, if you want all of
		them to be, for example, red, add this to the <style> area:
			circle{ fill: red; stroke: red }
		
		If you just want a certain group of dots coloured, you can group them
		by adding a "class"-attribute like this:
			<circle cx="500" cy="100" r="2" class="reddot" />
		And add the following rule to the <style> area:
			circle.reddot{ fill: red; stroke: red; }
		Only the circles who belong to the "reddot" class will turn red 

		you can use rgb values as well: fill: rgb(255,0,0);

*/
		fprintf(svgofp,"    <circle cx=\"%.10f\" cy=\"%.10f\" r=\"2\" ", 
			(0.5*w1 + w2)*1000, -(w1*866.0254038));
		fprintf(svgofp,"fill=\"rgb(%d%%, %d%%, %d%%)\" />\n", (int)(100*red), (int)(100*green), (int)(100*blue));
	}
} /* plotlmpointcolor */
#endif




/* last lines of EPSF likelihood mapping file */
//void finisheps(FILE *ofp, unsigned long **countarr)
void finisheps(FILE *ofp, vector<SeqQuartetInfo> lmap_seq_quartet_info, int leafNum, unsigned long Numquartets)
{
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "%% second triangle (the one with 3 basins)\n");
	fprintf(ofp, "/secondtriangle {\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.0 tl 0.0 tl moveto\n");
	fprintf(ofp, " 1.0 tl 0.0 tl lineto\n");
	fprintf(ofp, " 0.5 tl 0.8660254038 tl lineto\n");
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.50 tl 0.2886751346 tl moveto\n");
	fprintf(ofp, " 0.50 tl 0.0000000000 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.50 tl 0.2886751346 tl moveto\n");
	fprintf(ofp, " 0.25 tl 0.4330127019 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.50 tl 0.2886751346 tl moveto\n");
	fprintf(ofp, " 0.75 tl 0.4330127019 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "0.44 tl 0.5 tl moveto %% up\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_AR1]*100.0/Numquartets);
	fprintf(ofp, "0.25 tl 0.15 tl moveto %% down left\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_AR3]*100.0/Numquartets);
	fprintf(ofp, "0.63 tl 0.15 tl moveto %% down right\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_AR2]*100.0/Numquartets);
	fprintf(ofp, "} def\n");
	fprintf(ofp, "%% third triangle (the one with 7 basins)\n");
	fprintf(ofp, "/thirdtriangle {\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.0 tl 0.0 tl moveto\n");
	fprintf(ofp, " 1.0 tl 0.0 tl lineto\n");
	fprintf(ofp, " 0.5 tl 0.8660254038 tl lineto\n");
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.25 tl 0.1443375673 tl moveto\n");
	fprintf(ofp, " 0.75 tl 0.1443375673 tl lineto\n");
	fprintf(ofp, " 0.50 tl 0.5773502692 tl lineto\n");
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.125 tl 0.2165063509 tl moveto\n");
	fprintf(ofp, " 0.250 tl 0.1443375673 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.375 tl 0.6495190528 tl moveto\n");
	fprintf(ofp, " 0.500 tl 0.5773502692 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.625 tl 0.6495190528 tl moveto\n");
	fprintf(ofp, " 0.500 tl 0.5773502692 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.875 tl 0.2165063509 tl moveto\n");
	fprintf(ofp, " 0.750 tl 0.1443375673 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.750 tl 0.00 tl moveto\n");
	fprintf(ofp, " 0.750 tl 0.1443375673 tl lineto\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, " 0.250 tl 0.00 tl moveto\n");
	fprintf(ofp, " 0.250 tl 0.1443375673 tl lineto\n");
	fprintf(ofp, "stroke\n");
	/* resolved quartets, top */
	fprintf(ofp, "0.42 tl 0.66 tl moveto %% up\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_REG1]*100.0/Numquartets);
	/* resolved quartets, bottom left */
	fprintf(ofp, "0.07 tl 0.05 tl moveto %% down left\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_REG3]*100.0/Numquartets);
	/* resolved quartets, bottom right */
	fprintf(ofp, "0.77 tl 0.05 tl moveto %% down right\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_REG2]*100.0/Numquartets);
	/* partly resolved quartets, bottom */
	fprintf(ofp, "0.43 tl 0.05 tl moveto %% down side\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_REG5]*100.0/Numquartets);
	/* unresolved quartets */
	fprintf(ofp, "0.43 tl 0.28 tl moveto %% center\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_REG7]*100.0/Numquartets);
	/* partly resolved quartets, top right */
	fprintf(ofp, "gsave\n");
	fprintf(ofp, "-60 rotate\n");
	fprintf(ofp, "-0.07 tl 0.77 tl moveto %% right side\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_REG4]*100.0/Numquartets);
	fprintf(ofp, "grestore\n");
	/* partly resolved quartets, top left */
	fprintf(ofp, "gsave\n");
	fprintf(ofp, "60 rotate\n");
	fprintf(ofp, "0.4 tl -0.09 tl moveto %% left side\n");
	fprintf(ofp, "(%.1f%%) show\n", (double) lmap_seq_quartet_info[leafNum].countarr[LM_REG6]*100.0/Numquartets);
	fprintf(ofp, "grestore\n");
	fprintf(ofp, "} def\n");
	fprintf(ofp, "%% print the other two triangles\n");
	fprintf(ofp, "-0.6 tl -1.2 tl 0.8660254038 mul translate\n");
	fprintf(ofp, "secondtriangle\n");
	fprintf(ofp, "1.2 tl 0 translate\n");
	fprintf(ofp, "thirdtriangle\n");	
	fprintf(ofp, "showpage\n");
	fprintf(ofp, "%%%%EOF\n");
} /* finisheps */

/*  Likelihood mapping to EPS file  */



/* switching on different color for each simplex area */
#define COL_TEST
#undef COL_TEST

#if 0  /*** old code from TREE-PUZZLE ***/
/* computes LM point from the three log-likelihood values,
   plots the point, and does some statistics */
/* input: log-likelihoods=b1/b2/b3, internal quartet edge length=bl1/bl2/bl3 */
/*        min edgelen=minlen, max edgelen=maxlen, how to color=coltype */
void makelmpoint(
	FILE *fpeps, 				/* output file pointer */
	FILE *fpsvg, 				/* output file pointer */
	double  b1, 				/* log-likelihood: ab|cd */
	double  b2, 				/* log-likelihood: ac|bd */
	double  b3, 				/* log-likelihood: ad|bc */
	int t1, int t2, int t3, int t4,		/* four taxa in the quartet */
	vector<SeqQuartetInfo> lmap_seq_quartet_info,				/* taxon-specific  */
	// unsigned long **countarr,				/* taxon-specific  */
	double  bl1, 				/* internal branchlength: ab|cd */
	double  bl2, 				/* internal branchlength: ac|bd */
	double  bl3, 				/* internal branchlength: ad|bc */
	double  minlen, 			/* lower branch-length bound */
	double  maxlen, 			/* upper branch-length bound */
	int     coltype				/* color to plot points on the branch-length bound */
)
{
	double w1, w2, w3, temp;
	unsigned char qpbranching;
	double temp1, temp2, temp3, onethird;
	double templog;
	unsigned char discreteweight[3], treebits[3];
	int nummin, nummax;
	
	onethird = 1.0/3.0;
	treebits[0] = (unsigned char) 1;
	treebits[1] = (unsigned char) 2;
	treebits[2] = (unsigned char) 4;

	/* sort in descending order */
	qweight[0] = b1;
	qweight[1] = b2;
	qweight[2] = b3;
	sort3doubles(qweight, qworder);

	/* compute Bayesian weights */
	templog = qweight[qworder[1]]-qweight[qworder[0]];
	if(templog < -TP_MAX_EXP_DIFF)	/* possible, since 1.0+exp(>36) == 1.0 */
	   qweight[qworder[1]] = 0.0;
	else
	   qweight[qworder[1]] = exp(templog);

        templog = qweight[qworder[2]]-qweight[qworder[0]];
	if(templog < -TP_MAX_EXP_DIFF)	/* possible, since 1.0+exp(>36) == 1.0 */
	   qweight[qworder[2]] = 0.0;
	else
	   qweight[qworder[2]] = exp(templog);

	qweight[qworder[0]] = 1.0;

	temp = qweight[0] + qweight[1] + qweight[2];
	qweight[0] = qweight[0]/temp;
	qweight[1] = qweight[1]/temp;
	qweight[2] = qweight[2]/temp;

	/* plot one point in likelihood mapping triangle */
	w1 = qweight[0];
	w2 = qweight[1];
	w3 = qweight[2];

	/* check areas 1,2,3 */
	if (treebits[qworder[0]] == 1) {
		countarr[Maxspc][LM_AR1]++;
		countarr[t1][LM_AR1]++;
		countarr[t2][LM_AR1]++;
		countarr[t3][LM_AR1]++;
		countarr[t4][LM_AR1]++;
	} else {
		if (treebits[qworder[0]] == 2) {
			countarr[Maxspc][LM_AR2]++;
			countarr[t1][LM_AR2]++;
			countarr[t2][LM_AR2]++;
			countarr[t3][LM_AR2]++;
			countarr[t4][LM_AR2]++;
		} else {
			countarr[Maxspc][LM_AR3]++;
			countarr[t1][LM_AR3]++;
			countarr[t2][LM_AR3]++;
			countarr[t3][LM_AR3]++;
			countarr[t4][LM_AR3]++;
		}
	}

	/* check out regions 1,2,3,4,5,6,7 */

	/* 100 distribution */
	temp1 = 1.0 - qweight[qworder[0]];
	sqdiff[0] = temp1*temp1 +
		qweight[qworder[1]]*qweight[qworder[1]] +
		qweight[qworder[2]]*qweight[qworder[2]];
	discreteweight[0] = treebits[qworder[0]];

	/* 110 distribution */
	temp1 = 0.5 - qweight[qworder[0]];
	temp2 = 0.5 - qweight[qworder[1]];
	sqdiff[1] = temp1*temp1 + temp2*temp2 +
		qweight[qworder[2]]*qweight[qworder[2]];
	discreteweight[1] = treebits[qworder[0]] + treebits[qworder[1]];

	/* 111 distribution */
	temp1 = onethird - qweight[qworder[0]];
	temp2 = onethird - qweight[qworder[1]];
	temp3 = onethird - qweight[qworder[2]];
	sqdiff[2] = temp1 * temp1 + temp2 * temp2 + temp3 * temp3;
	discreteweight[2] = (unsigned char) 7;

	/* sort in descending order */
	sort3doubles(sqdiff, sqorder);
			
	qpbranching = (unsigned char) discreteweight[sqorder[2]];
							
/*
                   (a,b)-(c,d)   => numclust == 4
                   (a,b)-(c,c)   => numclust == 3
                   (a,a)-(b,b)   => numclust == 2

                     1l/\1r
                      /  \
                     /  1 \
                 6u / \  / \ 4u
                   /   \/   \
               6d /    /\    \ 4d
                 / 6  /  \  4 \
                /\   /  7 \   /\
             3u/  \ /______\ /  \2u
              / 3  |    5   |  2 \
             /_____|________|_____\
	        3d   5l  5r   2d
    (a,d)-(b,c)                  (a,c)-(b,d)   => numclust == 4
    (a,c)-(b,c)                  (a,c)-(b,c)   => numclust == 3
    (a,b)-(a,b)                  (a,b)-(a,b)   => numclust == 2

*/

	if (qpbranching == 1) {
		countarr[Maxspc][LM_REG1]++;
		countarr[t1][LM_REG1]++;
		countarr[t2][LM_REG1]++;
		countarr[t3][LM_REG1]++;
		countarr[t4][LM_REG1]++;
	}
	if (qpbranching == 2) {
		countarr[Maxspc][LM_REG2]++;
		countarr[t1][LM_REG2]++;
		countarr[t2][LM_REG2]++;
		countarr[t3][LM_REG2]++;
		countarr[t4][LM_REG2]++;
	}
	if (qpbranching == 4) {
		countarr[Maxspc][LM_REG3]++;
		countarr[t1][LM_REG3]++;
		countarr[t2][LM_REG3]++;
		countarr[t3][LM_REG3]++;
		countarr[t4][LM_REG3]++;
	}
	if (qpbranching == 3) {
		countarr[Maxspc][LM_REG4]++;
		countarr[t1][LM_REG4]++;
		countarr[t2][LM_REG4]++;
		countarr[t3][LM_REG4]++;
		countarr[t4][LM_REG4]++;
	}
	if (qpbranching == 6) {
		countarr[Maxspc][LM_REG5]++;
		countarr[t1][LM_REG5]++;
		countarr[t2][LM_REG5]++;
		countarr[t3][LM_REG5]++;
		countarr[t4][LM_REG5]++;
	}
	if (qpbranching == 5) {
		countarr[Maxspc][LM_REG6]++;
		countarr[t1][LM_REG6]++;
		countarr[t2][LM_REG6]++;
		countarr[t3][LM_REG6]++;
		countarr[t4][LM_REG6]++;
	}
	if (qpbranching == 7) {
		countarr[Maxspc][LM_REG7]++;
		countarr[t1][LM_REG7]++;
		countarr[t2][LM_REG7]++;
		countarr[t3][LM_REG7]++;
		countarr[t4][LM_REG7]++;
	}

#if 0
	switch (coltype) {
		case LMAP_AREA:
			switch (qpbranching) {
				case 1: /* 001 */
					plotlmpointcolor(fpeps, fpsvg, w1, w2, 0, 0, 1);
					break;
				case 2: /* 010 */
					plotlmpointcolor(fpeps, fpsvg, w1, w2, 0, 1, 0);
					break;
				case 4: /* 100 */
					plotlmpointcolor(fpeps, fpsvg, w1, w2, 1, 0, 0);
					break;
				case 3: /* 011 */
					plotlmpointcolor(fpeps, fpsvg, w1, w2, 0, 1, 1);
					break;
				case 5: /* 101 */
					plotlmpointcolor(fpeps, fpsvg, w1, w2, 1, 0, 1);
					break;
				case 6: /* 110 */
					plotlmpointcolor(fpeps, fpsvg, w1, w2, 1, 1, 0);
					break;
				case 7: /* 111 */
					plotlmpoint(fpeps, fpsvg, w1, w2);
					break;
			}
			break;
		case LMAP_MINMAXLEN:
			nummin = nummax = 0;
			if (bl1 <= minlen) nummin++;
			if (bl2 <= minlen) nummin++;
			if (bl3 <= minlen) nummin++;
			if (nummin == 2) {
				plotlmpointcolor(fpeps, fpsvg, w1, w2, 1, 0, 0);
			} else { 
				if (nummin == 1) {
					plotlmpointcolor(fpeps, fpsvg, w1, w2, 0, 0, 1);
				} else {
					plotlmpoint(fpeps, fpsvg, w1, w2);
				}
			}
			break;
		case LMAP_NORMAL:
		default:
#endif
			plotlmpoint(fpeps, fpsvg, w1, w2);
#if 0
			break;
	}
#endif
} /* makelmpoint */
#endif

/****************************************/
/*  end of Likelihood mapping routines  */
/****************************************/


/***************************************************************/ 


//*** end of likelihood mapping stuff (imported from TREE-PUZZLE's lmap.c) (HAS)


void PhyloTree::computeQuartetLikelihoods(vector<QuartetInfo> &lmap_quartet_info, QuartetGroups &LMGroups) {

    if (leafNum <= 4) 
        outError("Tree must have 5 or more taxa with unique sequences!");
        
    lmap_quartet_info.resize(params->lmap_num_quartets);
    
    int qc[] = {0, 1, 2, 3,  0, 2, 1, 3,  0, 3, 1, 2};
    
    double onethird = 1.0/3.0;
    unsigned char treebits[] = {1, 2, 4};

    int sizeA, sizeB, sizeC, sizeD, numGroups;
    int size3, size2, size1, size0;

LMGroups.numGroups = 0;
    if(LMGroups.numGroups == 0) { /* no grouping */
	LMGroups.numGroups = 1;
	LMGroups.GroupA.resize(leafNum);
    	for (int s = 0; s<leafNum; s++) LMGroups.GroupA[s] = s;
	LMGroups.numGrpSeqs[0] = leafNum; /* cluster A */
	LMGroups.numGrpSeqs[1] = 0; /* cluster B */
	LMGroups.numGrpSeqs[2] = 0; /* cluster C */
	LMGroups.numGrpSeqs[3] = 0; /* cluster D */
	LMGroups.numGrpSeqs[4] = 0; /* excluded */
	LMGroups.numQuartSeqs  = leafNum; /* all sequences in analysis */
	LMGroups.numSeqs       = leafNum; /* all sequences in alignment */
	LMGroups.numGroups = 1;
    }

    numGroups = LMGroups.numGroups;
    sizeA = LMGroups.numGrpSeqs[0]; /* cluster A */
    sizeB = LMGroups.numGrpSeqs[1]; /* cluster B */
    sizeC = LMGroups.numGrpSeqs[2]; /* cluster C */
    sizeD = LMGroups.numGrpSeqs[3]; /* cluster D */

    switch(LMGroups.numGroups){
	case 1: 
	   if(sizeA <= 4) 
		outError("Likelihood Mapping requires 5 or more taxa with unique sequences!"); 
	   break;
	case 2: 
	   if((sizeA < 2)||(sizeB < 2)) 
		outError("2-cluster Likelihood Mapping requires clusters A and B to have >=2 taxa with unique sequences!"); 
	   break;
	case 3: 
	   if((sizeA < 2)||(sizeB < 1)||(sizeC < 1)) 
		outError("3-cluster Likelihood Mapping requires cluster A to have >=2 and clusters B and C >=1 taxa with unique sequences!"); 
	   break;
	case 4: 
	   if((sizeA < 1)||(sizeB < 1)||(sizeA < 1)||(sizeD < 1)) 
		outError("4-cluster Likelihood Mapping requires all 4 clusters to have >0 taxa with unique sequences!"); 
	   break;
	default: 
	   outError("Unknown Likelihood Mapping mode! PLEASE report this to the developers!"); 
	   break;
    }
    
    switch(LMGroups.numGroups){
	case 1: 
	   size3 = sizeA-3;
	   size2 = sizeA-2;
	   size1 = sizeA-1;
	   size0 = sizeA;
	   LMGroups.uniqueQuarts = 1     + size3 +
	                           size2 * (size2-1) / 2 +
	                           size1 * (size1-1) * (size1-2) / 6 +
	                           size0 * (size0-1) * (size0-2) * (size0-3) / 24;
	   break;
	case 2: 
	   LMGroups.uniqueQuarts = (sizeA * (sizeA - 1)) / 2 * (sizeB * (sizeB - 1)) / 2; break;
	case 3: 
	   LMGroups.uniqueQuarts = (sizeA * (sizeA - 1)) / 2 * sizeB * sizeC; break;
	case 4: 
	   LMGroups.uniqueQuarts = sizeA * sizeB * sizeC * sizeD; break;
	default: 
	   outError("Unknown Likelihood Mapping mode! PLEASE report this to the developers!"); 
	   break;
    }
    

#ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
#endif
    for (int qid = 0; qid < params->lmap_num_quartets; qid++) { /*** draw lmap_num_quartets quartets randomly ***/
	// fprintf(stderr, "%d\n", qid); 

        // uniformly draw 4 taxa
	// (a) sample taxon 1
        // was: lmap_quartet_info[qid].seqID[0] = random_int(leafNum);
	lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[random_int(sizeA)];

        do {
	    // (b) sample taxon 2 according to the number of clusters
            // was: lmap_quartet_info[qid].seqID[1] = random_int(leafNum);
	    switch(numGroups){
		case 1: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[random_int(sizeA)]; break; // a1,A2|a3,a4
		case 2: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[random_int(sizeA)]; break; // a1,A2|b1,b2
		case 3: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[random_int(sizeA)]; break; // a1,A2|b, c
		case 4: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupB[random_int(sizeB)]; break; // a ,B |c, d
		default: outError("Unknown Likelihood Mapping sampling mode! PLEASE report this to the developers!"); break;
	    }
        } while (lmap_quartet_info[qid].seqID[1] == lmap_quartet_info[qid].seqID[0]);
        do {
	    // (c) sample taxon 3 according to the number of clusters
            // was: lmap_quartet_info[qid].seqID[2] = random_int(leafNum);
	    switch(numGroups){
		case 1: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupA[random_int(sizeA)]; break; // a1,a2|A3,a4
		case 2: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupB[random_int(sizeB)]; break; // a1,a2|B1,b2
		case 3: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupB[random_int(sizeB)]; break; // a1,a2|B, c
		case 4: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupC[random_int(sizeC)]; break; // a ,b |C, d
		default: outError("Unknown Likelihood Mapping sampling mode! PLEASE report this to the developers!"); break;
	    }
        } while (lmap_quartet_info[qid].seqID[2] == lmap_quartet_info[qid].seqID[0] || lmap_quartet_info[qid].seqID[2] == lmap_quartet_info[qid].seqID[1]);
        do {
	    // (d) sample taxon 4 according to the number of clusters
            // was: lmap_quartet_info[qid].seqID[3] = random_int(leafNum);
	    switch(numGroups){
		case 1: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupA[random_int(sizeA)]; break; // a1,a2|a3,A4
		case 2: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupB[random_int(sizeB)]; break; // a1,a2|b1,B2
		case 3: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupC[random_int(sizeC)]; break; // a1,a2|b, C
		case 4: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupD[random_int(sizeD)]; break; // a ,b |c, D
		default: outError("Unknown Likelihood Mapping sampling mode! PLEASE report this to the developers!"); break;
	    }
        } while (lmap_quartet_info[qid].seqID[3] == lmap_quartet_info[qid].seqID[0] || lmap_quartet_info[qid].seqID[3] == lmap_quartet_info[qid].seqID[1]
            || lmap_quartet_info[qid].seqID[3] == lmap_quartet_info[qid].seqID[2]);
            
// fprintf(stderr, "qqq%d: %d, %d, %d, %d\n", qid, lmap_quartet_info[qid].seqID[0], lmap_quartet_info[qid].seqID[1], lmap_quartet_info[qid].seqID[2], lmap_quartet_info[qid].seqID[3]); 

	// *** taxa should not be sorted, because that changes the corners a dot is assigned to - removed HAS ;^)
        // obsolete: sort(lmap_quartet_info[qid].seqID, lmap_quartet_info[qid].seqID+4); // why sort them?!? HAS ;^)
            
        // initialize sub-alignment and sub-tree
        Alignment *quartet_aln;
        PhyloTree *quartet_tree;
        if (aln->isSuperAlignment()) {
            quartet_aln = new SuperAlignment;
        } else {
            quartet_aln = new Alignment;
        }
        IntVector seq_id;
        seq_id.insert(seq_id.begin(), lmap_quartet_info[qid].seqID, lmap_quartet_info[qid].seqID+4);
        quartet_aln->extractSubAlignment(aln, seq_id, 0);
        if (isSuperTree()) {
            quartet_tree = new PhyloSuperTree((SuperAlignment*)quartet_aln, (PhyloSuperTree*)this);
        } else {
            quartet_tree = new PhyloTree(quartet_aln);
        }

        // set up parameters
        quartet_tree->setParams(params);
        quartet_tree->optimize_by_newton = params->optimize_by_newton;
        quartet_tree->setLikelihoodKernel(params->SSE);
        
        // set model and rate
        quartet_tree->setModelFactory(model_factory);
        quartet_tree->setModel(getModel());
        quartet_tree->setRate(getRate());
        // NOTE: we don't need to set phylo_tree in model and rate because parameters are not reoptimized
        
        
        
        // loop over 3 quartets to compute likelihood
        for (int k = 0; k < 3; k++) {
            string quartet_tree_str;
            quartet_tree_str = "(" + quartet_aln->getSeqName(qc[k*4]) + "," + quartet_aln->getSeqName(qc[k*4+1]) + ",(" + 
                quartet_aln->getSeqName(qc[k*4+2]) + "," + quartet_aln->getSeqName(qc[k*4+3]) + "));";
            quartet_tree->readTreeStringSeqName(quartet_tree_str);
            quartet_tree->initializeAllPartialLh();
            quartet_tree->wrapperFixNegativeBranch(true);
            // optimize branch lengths with logl_epsilon=0.1 accuracy
            lmap_quartet_info[qid].logl[k] = quartet_tree->optimizeAllBranches(10, 0.1);
        }
        // reset model & rate so that they are not deleted
        quartet_tree->setModel(NULL);
        quartet_tree->setModelFactory(NULL);
        quartet_tree->setRate(NULL);

        delete quartet_tree;
        delete quartet_aln;

	// determine likelihood order
	int qworder[3]; // local (thread-safe) vector for sorting

	if (lmap_quartet_info[qid].logl[0] > lmap_quartet_info[qid].logl[1]) {
		if(lmap_quartet_info[qid].logl[2] > lmap_quartet_info[qid].logl[0]) {
			qworder[0] = 2;
			qworder[1] = 0;
			qworder[2] = 1;		
		} else if (lmap_quartet_info[qid].logl[2] < lmap_quartet_info[qid].logl[1]) {
			qworder[0] = 0;
			qworder[1] = 1;
			qworder[2] = 2;		
		} else {
			qworder[0] = 0;
			qworder[1] = 2;
			qworder[2] = 1;		
		}
	} else {
		if(lmap_quartet_info[qid].logl[2] > lmap_quartet_info[qid].logl[1]) {
			qworder[0] = 2;
			qworder[1] = 1;
			qworder[2] = 0;		
		} else if (lmap_quartet_info[qid].logl[2] < lmap_quartet_info[qid].logl[0]) {
			qworder[0] = 1;
			qworder[1] = 0;
			qworder[2] = 2;		
		} else {
			qworder[0] = 1;
			qworder[1] = 2;
			qworder[2] = 0;		
		}
	}

	// compute Bayesian weights
	double temp;

	lmap_quartet_info[qid].qweight[0] = lmap_quartet_info[qid].logl[0];
	lmap_quartet_info[qid].qweight[1] = lmap_quartet_info[qid].logl[1];
	lmap_quartet_info[qid].qweight[2] = lmap_quartet_info[qid].logl[2];

	temp = lmap_quartet_info[qid].qweight[qworder[1]]-lmap_quartet_info[qid].qweight[qworder[0]];
	if(temp < -TP_MAX_EXP_DIFF)	/* possible, since 1.0+exp(>36) == 1.0 */
	   lmap_quartet_info[qid].qweight[qworder[1]] = 0.0;
	else
	   lmap_quartet_info[qid].qweight[qworder[1]] = exp(temp);

        temp = lmap_quartet_info[qid].qweight[qworder[2]]-lmap_quartet_info[qid].qweight[qworder[0]];
	if(temp < -TP_MAX_EXP_DIFF)	/* possible, since 1.0+exp(>36) == 1.0 */
	   lmap_quartet_info[qid].qweight[qworder[2]] = 0.0;
	else
	   lmap_quartet_info[qid].qweight[qworder[2]] = exp(temp);

	lmap_quartet_info[qid].qweight[qworder[0]] = 1.0;

	temp = lmap_quartet_info[qid].qweight[0] + lmap_quartet_info[qid].qweight[1] + lmap_quartet_info[qid].qweight[2];
	lmap_quartet_info[qid].qweight[0] = lmap_quartet_info[qid].qweight[0]/temp;
	lmap_quartet_info[qid].qweight[1] = lmap_quartet_info[qid].qweight[1]/temp;
	lmap_quartet_info[qid].qweight[2] = lmap_quartet_info[qid].qweight[2]/temp;

	// determine which of the three corners (only meaningful if seqIDs NOT sorted)
	if (treebits[qworder[0]] == 1) {
		lmap_quartet_info[qid].corner=0;
	} else {
		if (treebits[qworder[0]] == 2) {
			lmap_quartet_info[qid].corner=1;
		} else {
			lmap_quartet_info[qid].corner=2;
		}
	}

	// determine which of the 7 regions (only meaningful if seqIDs NOT sorted)
	double temp1, temp2, temp3;
	unsigned char discreteweight[3];
	double sqdiff[3];

	/* 100 distribution */
	temp1 = 1.0 - lmap_quartet_info[qid].qweight[qworder[0]];
	sqdiff[0] = temp1*temp1 +
		lmap_quartet_info[qid].qweight[qworder[1]]*lmap_quartet_info[qid].qweight[qworder[1]] +
		lmap_quartet_info[qid].qweight[qworder[2]]*lmap_quartet_info[qid].qweight[qworder[2]];
	discreteweight[0] = treebits[qworder[0]];

	/* 110 distribution */
	temp1 = 0.5 - lmap_quartet_info[qid].qweight[qworder[0]];
	temp2 = 0.5 - lmap_quartet_info[qid].qweight[qworder[1]];
	sqdiff[1] = temp1*temp1 + temp2*temp2 +
		lmap_quartet_info[qid].qweight[qworder[2]]*lmap_quartet_info[qid].qweight[qworder[2]];
	discreteweight[1] = treebits[qworder[0]] + treebits[qworder[1]];

	/* 111 distribution */
	temp1 = onethird - lmap_quartet_info[qid].qweight[qworder[0]];
	temp2 = onethird - lmap_quartet_info[qid].qweight[qworder[1]];
	temp3 = onethird - lmap_quartet_info[qid].qweight[qworder[2]];
	sqdiff[2] = temp1 * temp1 + temp2 * temp2 + temp3 * temp3;
	discreteweight[2] = (unsigned char) 7;

	/* sort in descending order */
	int sqorder[3]; // local (thread-safe) vector for sorting
	if (sqdiff[0] > sqdiff[1]) {
		if(sqdiff[2] > sqdiff[0]) {
			sqorder[0] = 2;
			sqorder[1] = 0;
			sqorder[2] = 1;		
		} else if (sqdiff[2] < sqdiff[1]) {
			sqorder[0] = 0;
			sqorder[1] = 1;
			sqorder[2] = 2;		
		} else {
			sqorder[0] = 0;
			sqorder[1] = 2;
			sqorder[2] = 1;		
		}
	} else {
		if(sqdiff[2] > sqdiff[1]) {
			sqorder[0] = 2;
			sqorder[1] = 1;
			sqorder[2] = 0;		
		} else if (sqdiff[2] < sqdiff[0]) {
			sqorder[0] = 1;
			sqorder[1] = 0;
			sqorder[2] = 2;		
		} else {
			sqorder[0] = 1;
			sqorder[1] = 2;
			sqorder[2] = 0;		
		}
	}


	// determine which of the 7 regions (only meaningful if seqIDs NOT sorted)
	unsigned char qpbranching = (unsigned char) discreteweight[sqorder[2]];

	if (qpbranching == 1) {
		lmap_quartet_info[qid].area=0; // LM_REG1 - top
	}
	if (qpbranching == 2) {
		lmap_quartet_info[qid].area=1; // LM_REG2 - right
	}
	if (qpbranching == 4) {
		lmap_quartet_info[qid].area=2; // LM_REG3 - left
	}

	if (qpbranching == 3) {
		lmap_quartet_info[qid].area=3; // LM_REG4
	}
	if (qpbranching == 6) {
		lmap_quartet_info[qid].area=4; // LM_REG5
	}
	if (qpbranching == 5) {
		lmap_quartet_info[qid].area=5; // LM_REG6
	}

	if (qpbranching == 7) {
		lmap_quartet_info[qid].area=6; // LM_REG7 - center 
	}

    } /*** end draw lmap_num_quartets quartets randomly ***/

} // end PhyloTree::computeQuartetLikelihoods


//**************************************

void PhyloTree::readLikelihoodMappingGroups() {
} // end PhyloTree::doLikelihoodMapping
//**************************************

void PhyloTree::doLikelihoodMapping() {
    // TODO For Heiko: Please add code here
    // vector<QuartetInfo> lmap_quartet_info;
    // vector<SeqQuartetInfo> lmap_seq_quartet_info;
    // int areacount[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // int cornercount[4] = {0, 0, 0, 0};
    int resolved, partly, unresolved;
    int qid;
    ofstream out;
    string filename;

    areacount[0] = 0;
    areacount[1] = 0;
    areacount[2] = 0;
    areacount[3] = 0;
    areacount[4] = 0;
    areacount[5] = 0;
    areacount[6] = 0;
    areacount[7] = 0;
    cornercount[0] = 0;
    cornercount[1] = 0;
    cornercount[2] = 0;
    cornercount[3] = 0;

    lmap_seq_quartet_info.resize(leafNum+1);
    for (qid = 0; qid < leafNum; qid++) {
        lmap_seq_quartet_info[qid].countarr[0] = 0;
        lmap_seq_quartet_info[qid].countarr[1] = 0;
        lmap_seq_quartet_info[qid].countarr[2] = 0;
        lmap_seq_quartet_info[qid].countarr[3] = 0;
        lmap_seq_quartet_info[qid].countarr[4] = 0;
        lmap_seq_quartet_info[qid].countarr[5] = 0;
        lmap_seq_quartet_info[qid].countarr[6] = 0;
        lmap_seq_quartet_info[qid].countarr[7] = 0;
        lmap_seq_quartet_info[qid].countarr[8] = 0;
        lmap_seq_quartet_info[qid].countarr[9] = 0;
    }

    computeQuartetLikelihoods(lmap_quartet_info, LMGroups);

    for (qid = 0; qid < params->lmap_num_quartets; qid++) {
	int tempreg;

	tempreg = lmap_quartet_info[qid].area;
        areacount[tempreg]++;
        lmap_seq_quartet_info[leafNum].countarr[tempreg]++; /* which of the 7 regions */
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[0]].countarr[tempreg]++;
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[1]].countarr[tempreg]++;
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[2]].countarr[tempreg]++;
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[3]].countarr[tempreg]++;

        tempreg = LM_AR1+lmap_quartet_info[qid].corner; /* which of the 3 corners */
        cornercount[lmap_quartet_info[qid].corner]++;
        lmap_seq_quartet_info[leafNum].countarr[tempreg]++;
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[0]].countarr[tempreg]++;
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[1]].countarr[tempreg]++;
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[2]].countarr[tempreg]++;
        lmap_seq_quartet_info[lmap_quartet_info[qid].seqID[3]].countarr[tempreg]++;
    }
    
    if (params->print_lmap_quartet_lh) {
        // print quartet file
        filename = (string)params->out_prefix + ".lmap.quartetlh";
        out.open(filename.c_str());
    }

    string lmap_svgfilename = (string)params->out_prefix + ".lmap.svg";
    FILE *svgout;
    svgout = fopen(lmap_svgfilename.c_str(), "w");
    initsvg(svgout);

    string lmap_epsfilename = (string)params->out_prefix + ".lmap.eps";
    FILE *epsout;
    epsout = fopen(lmap_epsfilename.c_str(), "w");
    initeps(epsout);

    for (qid = 0; qid < params->lmap_num_quartets; qid++) {

	plotlmpointeps(epsout, 
		// double w1, double w2)
                lmap_quartet_info[qid].qweight[0],
                lmap_quartet_info[qid].qweight[1]);

	plotlmpointsvg(svgout, 
		// double w1, double w2)
                lmap_quartet_info[qid].qweight[0],
                lmap_quartet_info[qid].qweight[1]);

    }

    if (params->print_lmap_quartet_lh) {
        // print quartet file
        for (qid = 0; qid < params->lmap_num_quartets; qid++) {
            out << "(" << lmap_quartet_info[qid].seqID[0] << ","
                << lmap_quartet_info[qid].seqID[1] << ","
                << lmap_quartet_info[qid].seqID[2] << ","
                << lmap_quartet_info[qid].seqID[3] << ")"
                << "\t" << lmap_quartet_info[qid].logl[0] 
                << "\t" << lmap_quartet_info[qid].logl[1] 
                << "\t" << lmap_quartet_info[qid].logl[2]
                << "\t" << lmap_quartet_info[qid].qweight[0] 
                << "\t" << lmap_quartet_info[qid].qweight[1] 
                << "\t" << lmap_quartet_info[qid].qweight[2] << endl;
        }

        PhyloTree::reportLikelihoodMapping(out);

    /**** begin of report output ****/
    /**** moved to PhyloTree::reportLikelihoodMapping ****/
#if 0
	// LM_REG1 0   /* top corner */
	// LM_REG2 1   /* bottom-right corner */
	// LM_REG3 2   /* bottom-left corner */
	// LM_REG4 3   /* right rectangle */
	// LM_REG5 4   /* bottom rectangle */
	// LM_REG6 5   /* left rectangle */
	// LM_REG7 6   /* center */
	// LM_AR1  7   /* top third */
	// LM_AR2  8   /* bottom-right third */
	// LM_AR3  9   /* bottom-left third */

	out << "LIKELIHOOD MAPPING ANALYSIS" << endl << endl;
	out << "Number of quartets: " << params->lmap_num_quartets << "(random choice)" << endl << endl;
	out << "Quartet trees are based on the selected model of substitution." << endl << endl;
	out << "Sequences are not grouped in clusters." << endl;

        out << endl << endl;
	out << "LIKELIHOOD MAPPING STATISTICS" << endl << endl;

	out << "           (a,b)-(c,d)                              (a,b)-(c,d)      " << endl;
	out << "                /\\                                      /\\           " << endl;
	out << "               /  \\                                    /  \\          " << endl;
	out << "              /    \\                                  /  1 \\         " << endl;
	out << "             /  a1  \\                                / \\  / \\        " << endl;
	out << "            /\\      /\\                              /   \\/   \\       " << endl;
	out << "           /  \\    /  \\                            /    /\\    \\      " << endl;
	out << "          /    \\  /    \\                          / 6  /  \\  4 \\     " << endl;
	out << "         /      \\/      \\                        /\\   /  7 \\   /\\    " << endl;
	out << "        /        |       \\                      /  \\ /______\\ /  \\   " << endl;
	out << "       /   a3    |    a2  \\                    / 3  |    5   |  2 \\  " << endl;
	out << "      /__________|_________\\                  /_____|________|_____\\ " << endl;
	out << "(a,d)-(b,c)            (a,c)-(b,d)      (a,b)-(c,d)            (a,c)-(b,d) " << endl << endl;
	out << "For more information about likelihood-mapping refer to" << endl;
   	out << "   Strimmer and von Haeseler (1997) PNAS 94:6815-6819" << endl;
	out << "   http://www.ncbi.nlm.nih.gov/pubmed/9192648" << endl;
	out << "and/or" << endl;
   	out << "   Schmidt and von Haeseler (2003) Current Protocols in Bioinformatics" << endl;
   	out << "   (by Baxevanis et al., Eds.), Unit 6, Wiley&Sons, New York." << endl;
	out << "   http://dx.doi.org/10.1002/0471250953.bi0606s17" << endl;


        out << endl << endl;
        out << "Quartet support of regions a1, a2, a3:" << endl << endl;
        out << "     #quartets    a1 (% a1)        a2 (% a2)        a3 (% a3)    name" << endl;
        for (qid = 0; qid < leafNum; qid++) {
	    //unsigned long sumq = lmap_seq_quartet_info[qid].countarr[0] + lmap_seq_quartet_info[qid].countarr[1] + lmap_seq_quartet_info[qid].countarr[2] + lmap_seq_quartet_info[qid].countarr[3] + lmap_seq_quartet_info[qid].countarr[4] + lmap_seq_quartet_info[qid].countarr[5] + lmap_seq_quartet_info[qid].countarr[6];
	    unsigned long sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];

            out.setf(ios::fixed, ios::floatfield); // set fixed floating format
            out.precision(2);
            out << setw(4) << qid+1 
                << setw(9) << sumq
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[7]
        	<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[7]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[8]
        	<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[8]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[9]
        	<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[9]/sumq << ")  " 
        	<< PhyloTree::aln->getSeqName(qid) << endl;
        }

        out << endl << endl << "Quartet support of areas 1-7:" << endl << endl;
        out << "                   resolved                                           partly                                             unresolved  name" << endl;
        out << "     #quartets     1 (% 1)          2 (% 2)          3 (% 3)          4 (% 4)          5 (% 5)          6 (% 6)          7 (% 7)" << endl;
        for (qid = 0; qid < leafNum; qid++) {
	    //unsigned long sumq = lmap_seq_quartet_info[qid].countarr[0] + lmap_seq_quartet_info[qid].countarr[1] + lmap_seq_quartet_info[qid].countarr[2] + lmap_seq_quartet_info[qid].countarr[3] + lmap_seq_quartet_info[qid].countarr[4] + lmap_seq_quartet_info[qid].countarr[5] + lmap_seq_quartet_info[qid].countarr[6];
	    unsigned long sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];

            out.setf(ios::fixed, ios::floatfield); // set fixed floating format
            out.precision(2);
            out << setw(4) << qid+1 
                << setw(9) << sumq
                << setw(7) << lmap_seq_quartet_info[qid].countarr[0] 
		<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[0]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[1] 
                << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[1]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[2]
                << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[2]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[3]
        	<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[3]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[4]
        	<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[4]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[5]
        	<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[5]/sumq << ") "
        	<< setw(7) << lmap_seq_quartet_info[qid].countarr[6]
        	<< " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[6]/sumq << ")  "
        	<< PhyloTree::aln->getSeqName(qid) << endl;
        }

        out << endl << endl << "Quartet resolution per sequence:" << endl << endl;
        out << "      #quartets   resolved           partly          unresolved  name" << endl;
        for (qid = 0; qid < leafNum; qid++) {
	    //unsigned long sumq = lmap_seq_quartet_info[qid].countarr[0] + lmap_seq_quartet_info[qid].countarr[1] + lmap_seq_quartet_info[qid].countarr[2] + lmap_seq_quartet_info[qid].countarr[3] + lmap_seq_quartet_info[qid].countarr[4] + lmap_seq_quartet_info[qid].countarr[5] + lmap_seq_quartet_info[qid].countarr[6];
	    unsigned long resolved = lmap_seq_quartet_info[qid].countarr[LM_REG1] + lmap_seq_quartet_info[qid].countarr[LM_REG2] + lmap_seq_quartet_info[qid].countarr[LM_REG3];
	    unsigned long partly   = lmap_seq_quartet_info[qid].countarr[LM_REG4] + lmap_seq_quartet_info[qid].countarr[LM_REG5] + lmap_seq_quartet_info[qid].countarr[LM_REG6];
	    unsigned long unres    = lmap_seq_quartet_info[qid].countarr[LM_REG7];
	    unsigned long sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];

            out.setf(ios::fixed, ios::floatfield); // set fixed floating format
            out.precision(2);
            out << setw(4) << qid+1 
                << setw(9) << sumq
                << setw(7) << resolved
		<< " (" << setw(6) << (double) 100.0*resolved/sumq << ") "
        	<< setw(7) << partly
                << " (" << setw(6) << (double) 100.0*partly/sumq << ") "
        	<< setw(7) << unres
                << " (" << setw(6) << (double) 100.0*unres/sumq << ")  "
        	<< PhyloTree::aln->getSeqName(qid) << endl;
        }
#endif
    }

    resolved   = areacount[0] + areacount[1] + areacount[2];
    partly     = areacount[3] + areacount[4] + areacount[5];
    unresolved = areacount[6];
	
    out << endl << "LIKELIHOOD MAPPING SUMMARY" << endl << endl;
    out << "Number of quartets: " << (resolved+partly+unresolved)
        << " (randomly drawn with replacement)" << endl << endl;
    out << "Overall quartet resolution:" << endl;
    out << "Number of fully resolved  quartets: " << resolved   
        << " (" << 100.0 * resolved/(resolved+partly+unresolved)   << "%)" << endl;
    out << "Number of partly resolved quartets: " << partly     
        << " (" << 100.0 * partly/(resolved+partly+unresolved)     << "%)" << endl;
    out << "Number of unresolved      quartets: " << unresolved 
        << " (" << 100.0 * unresolved/(resolved+partly+unresolved) << "%)" << endl << endl;

#if 0
#endif

    /**** end of report output ****/
    /**** moved to PhyloTree::reportLikelihoodMapping ****/


    if (params->print_lmap_quartet_lh) {
        // print quartet file
        out.close();        
        cout << "likelihood mapping results written to " << filename << endl;
    }

    finishsvg(svgout, lmap_seq_quartet_info, leafNum, params->lmap_num_quartets);
    fclose(svgout);        
    cout << "likelihood mapping plot (SVG) written to " << lmap_svgfilename << endl;

    finisheps(epsout, lmap_seq_quartet_info, leafNum, params->lmap_num_quartets);
    fclose(epsout);        
    cout << "likelihood mapping plot (EPS) written to " << lmap_epsfilename << endl;


//    cout << "\nOverall quartet resolution: (from " << (resolved+partly+unresolved) << " randomly drawn quartets)" << endl;
//    cout << "Fully resolved quartets:  " << resolved   << " (= "
//        << (double) resolved * 100.0   / (resolved+partly+unresolved) << "%)" << endl;
//    cout << "Partly resolved quartets: " << partly     << " (= "
//        << (double) partly * 100.0     / (resolved+partly+unresolved) << "%)" << endl;
//    cout << "Unresolved quartets:      " << unresolved << " (= "
//        << (double) unresolved * 100.0 / (resolved+partly+unresolved) << "%)" << endl << endl;

} // end PhyloTree::doLikelihoodMapping




void PhyloTree::reportLikelihoodMapping(ofstream &out) {
    // int areacount[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // int cornercount[4] = {0, 0, 0, 0};
    int resolved, partly, unresolved;
    int qid;
    int leafNum;
    leafNum = PhyloTree::aln->getNSeq();
    // vector<QuartetInfo> lmap_quartet_info;
    // vector<SeqQuartetInfo> lmap_seq_quartet_info;


	// LM_REG1 0   /* top corner */
	// LM_REG2 1   /* bottom-right corner */
	// LM_REG3 2   /* bottom-left corner */
	// LM_REG4 3   /* right rectangle */
	// LM_REG5 4   /* bottom rectangle */
	// LM_REG6 5   /* left rectangle */
	// LM_REG7 6   /* center */
	// LM_AR1  7   /* top third */
	// LM_AR2  8   /* bottom-right third */
	// LM_AR3  9   /* bottom-left third */
#if 0
#define LM_REG1 0   /* top corner */
#define LM_REG2 1   /* bottom-right corner */
#define LM_REG3 2   /* bottom-left corner */
#define LM_REG4 3   /* right rectangle */
#define LM_REG5 4   /* bottom rectangle */
#define LM_REG6 5   /* left rectangle */
#define LM_REG7 6   /* center */
#define LM_AR1  7   /* top third */
#define LM_AR2  8   /* bottom-right third */
#define LM_AR3  9   /* bottom-left third */
#endif

	out << "LIKELIHOOD MAPPING ANALYSIS" << endl;
	out << "---------------------------" << endl << endl;
	out << "Number of quartets: " << params->lmap_num_quartets << " (random choice)" << endl << endl;
	out << "Quartet trees are based on the selected model of substitution." << endl << endl;
	out << "Sequences are not grouped in clusters." << endl;

        out << endl << endl;
	out << "LIKELIHOOD MAPPING STATISTICS" << endl;
	out << "-----------------------------" << endl << endl;

	out << "           (a,b)-(c,d)                              (a,b)-(c,d)      " << endl;
	out << "                /\\                                      /\\           " << endl;
	out << "               /  \\                                    /  \\          " << endl;
	out << "              /    \\                                  /  1 \\         " << endl;
	out << "             /  a1  \\                                / \\  / \\        " << endl;
	out << "            /\\      /\\                              /   \\/   \\       " << endl;
	out << "           /  \\    /  \\                            /    /\\    \\      " << endl;
	out << "          /    \\  /    \\                          / 6  /  \\  4 \\     " << endl;
	out << "         /      \\/      \\                        /\\   /  7 \\   /\\    " << endl;
	out << "        /        |       \\                      /  \\ /______\\ /  \\   " << endl;
	out << "       /   a3    |    a2  \\                    / 3  |    5   |  2 \\  " << endl;
	out << "      /__________|_________\\                  /_____|________|_____\\ " << endl;
	out << "(a,d)-(b,c)            (a,c)-(b,d)      (a,b)-(c,d)            (a,c)-(b,d) " << endl << endl;
	out << "Division of the likelihood mapping plots into 3 or 7 areas." << endl;
	out << "On the left the areas show support for one of the different groupings like (a,b|c,d)." << endl;
	out << "On the right the right quartets falling into the areas 1, 2 and 3 are informative." << endl;
	out << "Those in the rectangles 4, 5 and 6 are partly informative and those in the center (7)" << endl;
	out << "are not informative." << endl;
	out << "Note, that the corners only make a difference if the sequences are clustered in groups." << endl << endl;
	out << "For more information about likelihood-mapping refer to" << endl;
   	out << " - Schmidt and von Haeseler (2009) The Phylogenetic Handbook, 2nd Ed." << endl;
   	out << "   (by Lemey et al., Eds.), 181-209, Cambridge Univ. Press, UK." << endl;
	out << "   http://www.thephylogenetichandbook.org" << endl;
   	out << " - Schmidt and von Haeseler (2003) Current Protocols in Bioinformatics" << endl;
   	out << "   (by Baxevanis et al., Eds.), Unit 6, Wiley&Sons, New York." << endl;
	out << "   http://dx.doi.org/10.1002/0471250953.bi0606s17" << endl;
	out << "and/or" << endl;
   	out << " - Strimmer and von Haeseler (1997) PNAS 94:6815-6819" << endl;
	out << "   http://www.ncbi.nlm.nih.gov/pubmed/9192648" << endl;


        out << endl << endl;
        out << "Quartet support of regions a1, a2, a3 (mainly for clustered analysis):" << endl << endl;
        out << "     #quartets    a1 (% a1)        a2 (% a2)        a3 (% a3)    name" << endl;
	out << "-----------------------------------------------------------------------------" << endl;
        for (qid = 0; qid <= leafNum; qid++) {
	    //unsigned long sumq = lmap_seq_quartet_info[qid].countarr[0] + lmap_seq_quartet_info[qid].countarr[1] + lmap_seq_quartet_info[qid].countarr[2] + lmap_seq_quartet_info[qid].countarr[3] + lmap_seq_quartet_info[qid].countarr[4] + lmap_seq_quartet_info[qid].countarr[5] + lmap_seq_quartet_info[qid].countarr[6];
	    unsigned long sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];

	    if (qid < leafNum) {
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << setw(4) << qid+1 
                   << setw(9) << sumq
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[7]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[7]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[8]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[8]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[9]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[9]/sumq << ")  " 
        	   << PhyloTree::aln->getSeqName(qid) << endl;
	    } else {
	       out << "-----------------------------------------------------------------------------" << endl;
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << "    "
                   << setw(9) << sumq
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[7]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[7]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[8]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[8]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[9]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[9]/sumq << ")  " << endl;
	    }
        }

        out << endl << endl << "Quartet support of areas 1-7 (mainly for clustered analysis):" << endl << endl;
        out << "                   resolved                                           partly                                             unresolved  name" << endl;
        out << "     #quartets     1 (% 1)          2 (% 2)          3 (% 3)          4 (% 4)          5 (% 5)          6 (% 6)          7 (% 7)" << endl;
	out << "------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
        for (qid = 0; qid <= leafNum; qid++) {
	    //unsigned long sumq = lmap_seq_quartet_info[qid].countarr[0] + lmap_seq_quartet_info[qid].countarr[1] + lmap_seq_quartet_info[qid].countarr[2] + lmap_seq_quartet_info[qid].countarr[3] + lmap_seq_quartet_info[qid].countarr[4] + lmap_seq_quartet_info[qid].countarr[5] + lmap_seq_quartet_info[qid].countarr[6];
	    unsigned long sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];

	    if (qid < leafNum) {
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << setw(4) << qid+1 
                   << setw(9) << sumq
                   << setw(7) << lmap_seq_quartet_info[qid].countarr[0] 
		   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[0]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[1] 
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[1]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[2]
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[2]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[3]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[3]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[4]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[4]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[5]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[5]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[6]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[6]/sumq << ")  "
        	   << PhyloTree::aln->getSeqName(qid) << endl;
	    } else {
	       out << "------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << "    "
                   << setw(9) << sumq
                   << setw(7) << lmap_seq_quartet_info[qid].countarr[0] 
		   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[0]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[1] 
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[1]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[2]
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[2]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[3]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[3]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[4]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[4]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[5]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[5]/sumq << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[6]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[6]/sumq << ")  " 
	           << endl << endl;
	    }
        }

        out << endl << endl << "Quartet resolution per sequence (phylogenetic information):" << endl << endl;
        out << "      #quartets   resolved           partly          unresolved  name" << endl;
	out << "-----------------------------------------------------------------------------" << endl;
        for (qid = 0; qid <= leafNum; qid++) {
	    //unsigned long sumq = lmap_seq_quartet_info[qid].countarr[0] + lmap_seq_quartet_info[qid].countarr[1] + lmap_seq_quartet_info[qid].countarr[2] + lmap_seq_quartet_info[qid].countarr[3] + lmap_seq_quartet_info[qid].countarr[4] + lmap_seq_quartet_info[qid].countarr[5] + lmap_seq_quartet_info[qid].countarr[6];
	    unsigned long resolved = lmap_seq_quartet_info[qid].countarr[LM_REG1] + lmap_seq_quartet_info[qid].countarr[LM_REG2] + lmap_seq_quartet_info[qid].countarr[LM_REG3];
	    unsigned long partly   = lmap_seq_quartet_info[qid].countarr[LM_REG4] + lmap_seq_quartet_info[qid].countarr[LM_REG5] + lmap_seq_quartet_info[qid].countarr[LM_REG6];
	    unsigned long unres    = lmap_seq_quartet_info[qid].countarr[LM_REG7];
	    unsigned long sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];

	    if (qid < leafNum) {
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << setw(4) << qid+1 
                   << setw(9) << sumq
                   << setw(7) << resolved
		   << " (" << setw(6) << (double) 100.0*resolved/sumq << ") "
        	   << setw(7) << partly
                   << " (" << setw(6) << (double) 100.0*partly/sumq << ") "
        	   << setw(7) << unres
                   << " (" << setw(6) << (double) 100.0*unres/sumq << ")  "
        	   << PhyloTree::aln->getSeqName(qid) << endl;
	    } else {
	       out << "-----------------------------------------------------------------------------" << endl;
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << "    "
                   << setw(9) << sumq
                   << setw(7) << resolved
		   << " (" << setw(6) << (double) 100.0*resolved/sumq << ") "
        	   << setw(7) << partly
                   << " (" << setw(6) << (double) 100.0*partly/sumq << ") "
        	   << setw(7) << unres
                   << " (" << setw(6) << (double) 100.0*unres/sumq << ")  " << endl;
	    }
        }

    resolved   = areacount[0] + areacount[1] + areacount[2];
    partly     = areacount[3] + areacount[4] + areacount[5];
    unresolved = areacount[6];
	
    // out << endl << "LIKELIHOOD MAPPING ANALYSIS" << endl << endl;
    // out << "Number of quartets: " << (resolved+partly+unresolved)
    //    << " (randomly drawn with replacement)" << endl << endl;
    out << "Overall quartet resolution:" << endl << endl;
    out << "Number of fully resolved  quartets (regions 1+2+3): " << resolved   
        << " (=" << 100.0 * resolved/(resolved+partly+unresolved)   << "%)" << endl;
    out << "Number of partly resolved quartets (regions 4+5+6): " << partly     
        << " (=" << 100.0 * partly/(resolved+partly+unresolved)     << "%)" << endl;
    out << "Number of unresolved      quartets (region 7)     : " << unresolved 
        << " (=" << 100.0 * unresolved/(resolved+partly+unresolved) << "%)" << endl << endl;

} // end PhyloTree::reportLikelihoodMapping
