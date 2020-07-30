//
//  quartet.cpp
//  iqtree
//
//  Created by Minh Bui on 24/07/15.
//
//

#include <stdio.h>
#include <string.h>

#include "phylotree.h"
#include "phylosupertree.h"
#include "model/partitionmodel.h"
#include "alignment/alignment.h"
#if 0 // (HAS-bla)
#include "tools.h"
#endif
#include "ncl/ncl.h"
#include "nclextra/msetsblock.h"
#include "nclextra/myreader.h"
// #include "lmap.c"

#ifdef _OPENMP
#include <omp.h>
#endif

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
void initsvg(FILE *ofp, QuartetGroups &LMGroups)
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
#endif /* LMAP_CLUSTER */
	if (LMGroups.numGroups == 2) { /* two cluster analysis */
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"500.0\"\n");
		fprintf(ofp,"	   y=\"-896.0254038\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_top_1\">(a,a)-(b,b)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_top_1\">(%s,%s)-(%s,%s)</text> <!-- (a,a|b,b) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[0]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[1]).c_str());
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"-30.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_left_1\">(a,b)-(a,b)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_left_1\">(%s,%s)-(%s,%s)</text> <!-- (a,b|a,b) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str());
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"1030.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_right_1\">(a,b)-(a,b)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_right_1\">(%s,%s)-(%s,%s)</text> <!-- (a,b|a,b) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str());
	}
	if (LMGroups.numGroups == 3) { /* three cluster analysis */
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"500.0\"\n");
		fprintf(ofp,"	   y=\"-896.0254038\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_top_1\">(a,b)-(c,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_top_1\">(%s,%s)-(%s,%s)</text> <!-- (a,b|c,c) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[2]).c_str(),(LMGroups.Name[2]).c_str());
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"-30.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_left_1\">(a,c)-(b,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_left_1\">(%s,%s)-(%s,%s)</text> <!-- (a,c|b,c) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[2]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[2]).c_str());
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"1030.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_right_1\">(a,c)-(b,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_right_1\">(%s,%s)-(%s,%s)</text> <!-- (a,c|b,c) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[2]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[2]).c_str());
	}
	if (LMGroups.numGroups == 4) { /* four cluster analysis */
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"500.0\"\n");
		fprintf(ofp,"	   y=\"-896.0254038\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_top_1\">(a,b)-(c,d)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_top_1\">(%s,%s)-(%s,%s)</text> <!-- (a,b|c,d) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[2]).c_str(),(LMGroups.Name[3]).c_str());
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"-30.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_left_1\">(a,d)-(b,c)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_left_1\">(%s,%s)-(%s,%s)</text> <!-- (a,d|b,c) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[3]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[2]).c_str());
		fprintf(ofp,"	<text\n");
		fprintf(ofp,"	   x=\"1030.0\"\n");
		fprintf(ofp,"	   y=\"60.0\"\n");
		fprintf(ofp,"	   text-anchor=\"middle\"\n");
		// fprintf(ofp,"	   id=\"label_right_1\">(a,c)-(b,d)</text> <!-- CHANGE HERE IF NECESSARY -->\n");
		fprintf(ofp,"	   id=\"label_right_1\">(%s,%s)-(%s,%s)</text> <!-- (a,c|b,d) - CHANGE HERE IF NECESSARY -->\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[2]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[3]).c_str());
	}

} /* initsvg */




void plotlmpointsvg(FILE *ofp, double w1, double w2)
{
	/* plot dots into triangle 1 (top) */
	fprintf(ofp,"	<circle cx=\"%.10f\" cy=\"%.10f\" r=\"2\" />\n", (0.5*w1 + w2)*1000, -(w1*866.0254038));
} /* plotlmpointsvg */



// void finishsvg(FILE *ofp, unsigned long **countarr)
void finishsvg(FILE *ofp, vector<SeqQuartetInfo> lmap_seq_quartet_info, int leafNum, int64_t Numquartets)
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
void initeps(FILE *ofp, QuartetGroups &LMGroups)
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
#endif /* LMAP_CLUSTER */
	if (LMGroups.numGroups == 2) { /* two cluster analysis */
		fprintf(ofp, "%% label corners\n");
		fprintf(ofp, "0.5 tl 0.9 tl moveto\n"); /* old: 0.375 0.9 */
		// fprintf(ofp, "((a,a)-(b,b)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,a|b,b) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[0]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[1]).c_str());
		fprintf(ofp, "-0.045 tl -0.08 tl moveto\n"); /* old: -0.16 -0.08 */
		// fprintf(ofp, "((a,b)-(a,b)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,b|a,b) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str());
		fprintf(ofp, "1.045 tl -0.08 tl moveto\n"); /* old: -0.92 -0.08 */
		// fprintf(ofp, "((a,b)-(a,b)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,b|a,b) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str());
	}
	if (LMGroups.numGroups == 3) { /* three cluster analysis */
		fprintf(ofp, "%% label corners\n");
		fprintf(ofp, "0.5 tl 0.9 tl moveto\n"); /* old: 0.375 0.9 */
		// fprintf(ofp, "((a,b)-(c,c)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,b|c,c) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[2]).c_str(),(LMGroups.Name[2]).c_str());
		fprintf(ofp, "-0.045 tl -0.08 tl moveto\n"); /* old: -0.16 -0.08 */
		// fprintf(ofp, "((a,c)-(b,c)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,c|b,c) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[2]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[2]).c_str());
		fprintf(ofp, "1.045 tl -0.08 tl moveto\n"); /* old: -0.92 -0.08 */
		// fprintf(ofp, "((a,c)-(b,c)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,c|b,c) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[2]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[2]).c_str());
	}
	if (LMGroups.numGroups == 4) { /* four cluster analysis */
		fprintf(ofp, "%% label corners\n");
		fprintf(ofp, "0.5 tl 0.9 tl moveto\n"); /* old: 0.375 0.9 */
		// fprintf(ofp, "((a,b)-(c,d)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,b|c,d) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[1]).c_str(),
			(LMGroups.Name[2]).c_str(),(LMGroups.Name[3]).c_str());
		fprintf(ofp, "-0.045 tl -0.08 tl moveto\n"); /* old: -0.16 -0.08 */
		// fprintf(ofp, "((a,d)-(b,c)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,d|b,c) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[3]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[2]).c_str());
		fprintf(ofp, "1.045 tl -0.08 tl moveto\n"); /* old: -0.92 -0.08 */
		// fprintf(ofp, "((a,c)-(b,d)) centershow %% CHANGE HERE IF NECESSARY\n");
		fprintf(ofp, "((%s,%s)-(%s,%s)) centershow %% (a,c|b,d) - CHANGE HERE IF NECESSARY\n",
			(LMGroups.Name[0]).c_str(),(LMGroups.Name[2]).c_str(),
			(LMGroups.Name[1]).c_str(),(LMGroups.Name[3]).c_str());
	}

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
void finisheps(FILE *ofp, vector<SeqQuartetInfo> lmap_seq_quartet_info, int leafNum, int64_t Numquartets)
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


/****************************************/
/*  end of Likelihood mapping routines  */
/****************************************/


/***************************************************************/ 


//*** end of likelihood mapping stuff (imported from TREE-PUZZLE's lmap.c) (HAS)


void PhyloTree::computeQuartetLikelihoods(vector<QuartetInfo> &lmap_quartet_info, QuartetGroups &LMGroups) {

    if (leafNum < 4) 
        outError("Tree must have 4 or more taxa with unique sequences!");
        
    int qc[] = {0, 1, 2, 3,  0, 2, 1, 3,  0, 3, 1, 2};
    
    double onethird = 1.0/3.0;
    unsigned char treebits[] = {1, 2, 4};

    int sizeA, sizeB, sizeC, sizeD, numGroups;
    int size3, size2, size1, size0;

	// LMGroups.numGroups = 0;
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
	   if(sizeA < 4) 
		outError("Likelihood Mapping requires 4 or more taxa with unique sequences!"); 
	   break;
	case 2: 
	   if((sizeA < 2)||(sizeB < 2)) 
		outError("2-cluster Likelihood Mapping requires clusters A and B to have >=2 taxa with unique sequences!"); 
	   break;
	case 3: 
	   if((sizeA < 1)||(sizeB < 1)||(sizeC < 2)) 
		outError("3-cluster Likelihood Mapping requires clusters B and C to have >=1 and cluster C >=2 taxa with unique sequences!"); 
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
	   size3 = sizeA-4;
	   size2 = sizeA-3;
	   size1 = sizeA-2;
	   size0 = sizeA-1;
	   LMGroups.uniqueQuarts = (int64_t)1 + size3 +
	                           (int64_t)size2 * (size2-1) / 2 +
	                           (int64_t)size1 * (size1-1) * (size1-2) / 6 +
	                           (int64_t)size0 * (size0-1) * (size0-2) * (size0-3) / 24;
	   break;
	case 2: 
	   LMGroups.uniqueQuarts = ((int64_t)sizeA * (sizeA - 1)) / 2 * (sizeB * (sizeB - 1)) / 2; break;
	case 3: 
	   LMGroups.uniqueQuarts = (int64_t)sizeA * sizeB * (sizeC * (sizeC - 1)) / 2; break;
	case 4: 
	   LMGroups.uniqueQuarts = (int64_t)sizeA * sizeB * sizeC * sizeD; break;
	default: 
	   outError("Unknown Likelihood Mapping mode! PLEASE report this to the developers!"); 
	   break;
    }

    if (params->lmap_num_quartets == 0)
        params->lmap_num_quartets = LMGroups.uniqueQuarts;
    if (params->lmap_num_quartets > LMGroups.uniqueQuarts) {
        cout << "INFO: Number of quartets is reduced to all unique quartets " << LMGroups.uniqueQuarts << endl; 
    }

    cout << "Computing " << params->lmap_num_quartets << " quartet likelihoods (one dot represents 100 quartets)." << endl << endl;
    
    lmap_quartet_info.resize(params->lmap_num_quartets);
    
    bool quartets_drawn = false;
    
    if (params->lmap_num_quartets == LMGroups.uniqueQuarts) {
        // draw all unique quartets now
        quartets_drawn = true;
        int64_t qid = 0;
        switch (numGroups) {
        case 1:
            for (int i0 = 0; i0 < sizeA-3; i0++)
                for (int i1 = i0+1; i1 < sizeA-2; i1++)
                    for (int i2 = i1+1; i2 < sizeA-1; i2++)
                        for (int i3 = i2+1; i3 < sizeA; i3++) {
                            lmap_quartet_info[qid].seqID[0] = LMGroups.GroupA[i0];
                            lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[i1];
                            lmap_quartet_info[qid].seqID[2] = LMGroups.GroupA[i2];
                            lmap_quartet_info[qid].seqID[3] = LMGroups.GroupA[i3];
                            qid++;
                        }
            break;

        case 2:
            for (int i0 = 0; i0 < sizeA-1; i0++)
                for (int i1 = i0+1; i1 < sizeA; i1++)
                    for (int i2 = 0; i2 < sizeB-1; i2++)
                        for (int i3 = i2+1; i3 < sizeB; i3++) {
                            lmap_quartet_info[qid].seqID[0] = LMGroups.GroupA[i0];
                            lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[i1];
                            lmap_quartet_info[qid].seqID[2] = LMGroups.GroupB[i2];
                            lmap_quartet_info[qid].seqID[3] = LMGroups.GroupB[i3];
                            qid++;
                        }
            break;

        case 3:
            for (int i0 = 0; i0 < sizeA; i0++)
                for (int i1 = 0; i1 < sizeB; i1++)
                    for (int i2 = 0; i2 < sizeC-1; i2++)
                        for (int i3 = i2+1; i3 < sizeC; i3++) {
                            lmap_quartet_info[qid].seqID[0] = LMGroups.GroupA[i0];
                            lmap_quartet_info[qid].seqID[1] = LMGroups.GroupB[i1];
                            lmap_quartet_info[qid].seqID[2] = LMGroups.GroupC[i2];
                            lmap_quartet_info[qid].seqID[3] = LMGroups.GroupC[i3];
                            qid++;
                        }
            break;
        case 4:
            for (int i0 = 0; i0 < sizeA; i0++)
                for (int i1 = 0; i1 < sizeB; i1++)
                    for (int i2 = 0; i2 < sizeC; i2++)
                        for (int i3 = 0; i3 < sizeD; i3++) {
                            lmap_quartet_info[qid].seqID[0] = LMGroups.GroupA[i0];
                            lmap_quartet_info[qid].seqID[1] = LMGroups.GroupB[i1];
                            lmap_quartet_info[qid].seqID[2] = LMGroups.GroupC[i2];
                            lmap_quartet_info[qid].seqID[3] = LMGroups.GroupD[i3];
                            qid++;
                        }
            break;

        default:
            break;
        }
        // sanity check
        ASSERT(qid == LMGroups.uniqueQuarts);
    }
    
    // fprintf(stderr,"XXX - #quarts: %d; #groups: %d, A: %d, B:%d, C:%d, D:%d\n", LMGroups.uniqueQuarts, LMGroups.numGroups, sizeA, sizeB, sizeC, sizeD);
    

#ifdef _OPENMP
    #pragma omp parallel
    {
    int *rstream;
    init_random(params->ran_seed + omp_get_thread_num(), false, &rstream);
#else
    int *rstream = randstream;
#endif    

#ifdef _OPENMP
    #pragma omp for schedule(guided)
#endif
    for (int64_t qid = 0; qid < params->lmap_num_quartets; qid++) { /*** draw lmap_num_quartets quartets randomly ***/
	// fprintf(stderr, "%I64d\n", qid); 

        // uniformly draw 4 taxa
	// (a) sample taxon 1
        // was: lmap_quartet_info[qid].seqID[0] = random_int(leafNum);
        if (!quartets_drawn) {
            // draw a random quartet
            lmap_quartet_info[qid].seqID[0] = LMGroups.GroupA[random_int(sizeA, rstream)];

            do {
            // (b) sample taxon 2 according to the number of clusters
                // was: lmap_quartet_info[qid].seqID[1] = random_int(leafNum);
            switch(numGroups){
            case 1: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[random_int(sizeA, rstream)]; break; // a1,A2|a3,a4
            case 2: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupA[random_int(sizeA, rstream)]; break; // a1,A2|b1,b2
            case 3: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupB[random_int(sizeB, rstream)]; break; // a ,B |c1,c2
            case 4: lmap_quartet_info[qid].seqID[1] = LMGroups.GroupB[random_int(sizeB, rstream)]; break; // a ,B |c, d
            default: outError("Unknown Likelihood Mapping sampling mode! PLEASE report this to the developers!"); break;
            }
            } while (lmap_quartet_info[qid].seqID[1] == lmap_quartet_info[qid].seqID[0]);
            do {
            // (c) sample taxon 3 according to the number of clusters
                // was: lmap_quartet_info[qid].seqID[2] = random_int(leafNum);
            switch(numGroups){
            case 1: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupA[random_int(sizeA, rstream)]; break; // a1,a2|A3,a4
            case 2: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupB[random_int(sizeB, rstream)]; break; // a1,a2|B1,b2
            case 3: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupC[random_int(sizeC, rstream)]; break; // a ,b |C1,c2
            case 4: lmap_quartet_info[qid].seqID[2] = LMGroups.GroupC[random_int(sizeC, rstream)]; break; // a ,b |C, d
            default: outError("Unknown Likelihood Mapping sampling mode! PLEASE report this to the developers!"); break;
            }
            } while (lmap_quartet_info[qid].seqID[2] == lmap_quartet_info[qid].seqID[0] || lmap_quartet_info[qid].seqID[2] == lmap_quartet_info[qid].seqID[1]);
            do {
            // (d) sample taxon 4 according to the number of clusters
                // was: lmap_quartet_info[qid].seqID[3] = random_int(leafNum);
            switch(numGroups){
            case 1: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupA[random_int(sizeA, rstream)]; break; // a1,a2|a3,A4
            case 2: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupB[random_int(sizeB, rstream)]; break; // a1,a2|b1,B2
            case 3: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupC[random_int(sizeC, rstream)]; break; // a ,b |c1,C2
            case 4: lmap_quartet_info[qid].seqID[3] = LMGroups.GroupD[random_int(sizeD, rstream)]; break; // a ,b |c, D
            default: outError("Unknown Likelihood Mapping sampling mode! PLEASE report this to the developers!"); break;
            }
            } while (lmap_quartet_info[qid].seqID[3] == lmap_quartet_info[qid].seqID[0] || lmap_quartet_info[qid].seqID[3] == lmap_quartet_info[qid].seqID[1]
                || lmap_quartet_info[qid].seqID[3] == lmap_quartet_info[qid].seqID[2]);
        }
// fprintf(stderr, "qqq%d: %d, %d, %d, %d\n", qid, lmap_quartet_info[qid].seqID[0], lmap_quartet_info[qid].seqID[1], lmap_quartet_info[qid].seqID[2], lmap_quartet_info[qid].seqID[3]); 

	// *** taxa should not be sorted, because that changes the corners a dot is assigned to - removed HAS ;^)
        // obsolete: sort(lmap_quartet_info[qid].seqID, lmap_quartet_info[qid].seqID+4); // why sort them?!? HAS ;^)

        // initialize sub-alignment and sub-tree
        Alignment *quartet_aln;
        if (aln->isSuperAlignment()) {
            quartet_aln = new SuperAlignment;
        } else {
            quartet_aln = new Alignment;
        }
        IntVector seq_id;
        seq_id.insert(seq_id.begin(), lmap_quartet_info[qid].seqID, lmap_quartet_info[qid].seqID+4);
        IntVector kept_partitions;
        // only keep partitions with at least 3 sequences
        quartet_aln->extractSubAlignment(aln, seq_id, 0, 3, &kept_partitions);
                
        if (kept_partitions.size() == 0) {
            // nothing kept
            for (int k = 0; k < 3; k++) {
                lmap_quartet_info[qid].logl[k] = -1.0;
            }
        } else {
            // something partition kept, do computations
            if (quartet_aln->ordered_pattern.empty())
                quartet_aln->orderPatternByNumChars(PAT_VARIANT);
            PhyloTree *quartet_tree;
            if (isSuperTree()) {
                quartet_tree = new PhyloSuperTree((SuperAlignment*)quartet_aln, (PhyloSuperTree*)this);
            } else {
                quartet_tree = new PhyloTree(quartet_aln);
            }

            // set up parameters
            quartet_tree->setParams(params);
            quartet_tree->optimize_by_newton = params->optimize_by_newton;
            quartet_tree->setLikelihoodKernel(params->SSE);
            quartet_tree->setNumThreads(num_threads);

            // set model and rate
            quartet_tree->setModelFactory(model_factory);
            quartet_tree->setModel(getModel());
            quartet_tree->setRate(getRate());

            // set up partition model
            if (isSuperTree()) {
                PhyloSuperTree *quartet_super_tree = (PhyloSuperTree*)quartet_tree;
                PhyloSuperTree *super_tree = (PhyloSuperTree*)this;
                for (int i = 0; i < quartet_super_tree->size(); i++) {
                    quartet_super_tree->at(i)->setModelFactory(super_tree->at(kept_partitions[i])->getModelFactory());
                    quartet_super_tree->at(i)->setModel(super_tree->at(kept_partitions[i])->getModel());
                    quartet_super_tree->at(i)->setRate(super_tree->at(kept_partitions[i])->getRate());
                    //quartet_super_tree->at(i)->aln->buildSeqStates(quartet_super_tree->at(i)->getModel()->seq_states);
                }
            } else {
                //quartet_aln->buildSeqStates(getModel()->seq_states);
            }
            
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

            if (isSuperTree()) {
                PhyloSuperTree *quartet_super_tree = (PhyloSuperTree*)quartet_tree;
                for (int i = 0; i < quartet_super_tree->size(); i++) {
                    quartet_super_tree->at(i)->setModelFactory(NULL);
                    quartet_super_tree->at(i)->setModel(NULL);
                    quartet_super_tree->at(i)->setRate(NULL);
                }
            }
            delete quartet_tree;
        }
        
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

	{
		int64_t count = (qid+1);
		if ((count % 100) == 0) {
			cout << '.';
			if ((count % 1000) == 0) { // separator after 10 dots
				cout << ' ';
				if ((count % 5000) == 0) // new-line after 50 dots
					cout << " : " << count << endl;
			}
			cout.flush();
		}
	}
    } /*** end draw lmap_num_quartets quartets randomly ***/
#ifdef _OPENMP
    finish_random(rstream);
    }
#endif

    if ((params->lmap_num_quartets % 5000) != 0) {
	cout << ". : " << params->lmap_num_quartets << flush << endl << endl;
    } else cout << endl;


    // restore seq_states
    /*
    if (isSuperTree()) {
        PhyloSuperTree *super_tree = (PhyloSuperTree*)this;
        for (int i = 0; i < super_tree->size(); i++) {
            super_tree->at(i)->aln->buildSeqStates(super_tree->at(i)->getModel()->seq_states);
        }
    } else {
        aln->buildSeqStates(getModel()->seq_states);
    }
     */
} // end PhyloTree::computeQuartetLikelihoods


//**************************************

/**
    read groups in following format "(A, B, C, D), (E, F), (G, H), (I);"
**/
void readGroupNewick(char *filename, MSetsBlock *sets_block) {
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename);
        in.exceptions(ios::badbit);
        char ch;
        string name;
        TaxaSetNameVector *sets = sets_block->getSets();
        while (!in.eof()) {
            // read the cluster
            TaxaSetName *myset = new TaxaSetName;
			sets->push_back(myset);            
            in >> ch;
            if (ch != '(')
                throw "Cluster does not start with '('";
            // read taxon name
            do {
                in >> ch;
                name = "";
                do {
                    name += ch;
                    ch = in.get();
                    if (controlchar(ch)) {
                        in >> ch;
                        break;
                    }
                } while (ch != ',' && ch != ')');
                renameString(name);
                myset->taxlist.push_back(name);
                if (ch == ',') continue; // continue to read next taxon name
                if (ch == ')') break;
                throw "No ',' or ')' found after " + name;
            } while (ch != ')');
            in >> ch;
            if (isalnum(ch)) {
                // read cluster name
                name = "";
                do {
                    name += ch;
                    ch = in.get();
                    if (controlchar(ch)) {
                        in >> ch;
                        break;
                    }
                } while (ch != ',' && ch != ';');
                myset->name = name;
            } else {
                myset->name = "Cluster" + convertIntToString(sets->size());
            }
            // check for duplicated name
            for (TaxaSetNameVector::iterator it = sets->begin(); it != sets->end()-1; it++)
                if ((*it)->name == myset->name)
                    throw "Duplicated cluster name " + myset->name;
            if (ch == ',') continue; // continue to read next cluster
            if (ch == ';') break;
            throw "No ',' or ';' found after " + name + ")";
        }
        
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (const char* str) {
        outError(str);
    } catch (string &str) {
        outError(str);
    }
}

void PhyloTree::readLikelihoodMappingGroups(char *filename, QuartetGroups &LMGroups) {

    int numsets, numtax, taxid;
    char clustchar;
    MSetsBlock *lmclusters;
    lmclusters = new MSetsBlock();
    cout << endl << "Reading likelihood mapping cluster file " << filename << "..." << endl;

    InputType intype = detectInputFile(filename);
    
    if (intype == IN_NEXUS) {
        MyReader nexus(filename);
        nexus.Add(lmclusters);
        MyToken token(nexus.inf);
        nexus.Execute(token);
    } else if (intype == IN_NEWICK) {
        readGroupNewick(filename, lmclusters);
    } else
        outError("Unknown cluster file format, please use NEXUS or RAxML-style format");

    if (lmclusters->getNSets() == 0)
        outError("No cluster found");

    // lmclusters->Report(cout);

    cout << "(The leading numbers represent the order from the master alignment.)" << endl << endl;

    TaxaSetNameVector *allsets = lmclusters->getSets();
    numsets = allsets->size();

    if(numsets > 5) outError("Only up to 4 Likelihood Mapping clusters allowed, plus one 'ignored' cluster!");

    int n = 0;
    for (TaxaSetNameVector::iterator i = allsets->begin(); i != allsets->end(); i++) {
	if ((*i)->name.compare("ignored")==0 || (*i)->name.compare("IGNORED")==0) {
		LMGroups.Name[4] = (*i)->name;
		numtax = (*i)->taxlist.size();
		LMGroups.numGrpSeqs[4] = numtax;
		LMGroups.GroupX.resize(numtax);
		cout << "Cluster \"" << LMGroups.Name[4] << "\" lists " << (*i)->taxlist.size() << " sequences to be ignored:" << endl;
		int t = 0;
		for (vector<string>::iterator it = (*i)->taxlist.begin(); it != (*i)->taxlist.end(); it++) {
			taxid = aln->getSeqID(*it);
			if (taxid < 0) {
				cout << "WARNING: unknown sequence name \"" << (*it) << "\"! Will be ignored." << endl;
			} else {
				LMGroups.GroupX[t] = taxid;
				// cout << "  " << (*it) << " (" << taxid << "," << LMGroups.GroupX[t] << ")" << endl;
				cout << "  " << LMGroups.GroupX[t]+1 << ". " << (*it) << endl;
				t++;
			}
		}
		if (numtax != t) {
			cout << "WARNING: ignored cluster did contain unknown sequence names!" << endl;
			LMGroups.numGrpSeqs[4] = t;
		}
	} else {
		switch(n){
			case 0: clustchar='A'; break;
			case 1: clustchar='B'; break;
			case 2: clustchar='C'; break;
			case 3: clustchar='D'; break;
    			default: outError("Only up to 4 Likelihood Mapping clusters allowed, plus one 'ignored' cluster!"); break;
		}
		LMGroups.Name[n] = (*i)->name;
		numtax = (*i)->taxlist.size();
		LMGroups.numGrpSeqs[n] = numtax;
		switch(n){
			case 0: LMGroups.GroupA.resize(numtax); break;
			case 1: LMGroups.GroupB.resize(numtax); break;
			case 2: LMGroups.GroupC.resize(numtax); break;
			case 3: LMGroups.GroupD.resize(numtax); break;
    			default: outError("Only up to 4 Likelihood Mapping clusters allowed, plus one 'ignored' cluster!"); break;
		}

		cout << "Cluster " << n+1 << " \"" << LMGroups.Name[n] << "\" lists " << (*i)->taxlist.size() << " sequences: " << endl;

		int t = 0;
		for (vector<string>::iterator it = (*i)->taxlist.begin(); it != (*i)->taxlist.end(); it++) {
			taxid = aln->getSeqID(*it);
			if (taxid < 0) {
				cout << "WARNING: sequence name \"" << (*it) << "\"! Will be ignored." << endl;
			} else {
				switch(n){
					case 0: LMGroups.GroupA[t] = taxid;
						cout << "  " << LMGroups.GroupA[t]+1 << ". " << (*it) << endl;
						break;
					case 1: LMGroups.GroupB[t] = taxid;
						cout << "  " << LMGroups.GroupB[t]+1 << ". " << (*it) << endl;
						break;
					case 2: LMGroups.GroupC[t] = taxid;
						cout << "  " << LMGroups.GroupC[t]+1 << ". " << (*it) << endl;
						break;
					case 3: LMGroups.GroupD[t] = taxid;
						cout << "  " << LMGroups.GroupD[t]+1 << ". " << (*it) << endl;
						break;
    					default: outError("Only up to 4 Likelihood Mapping clusters allowed, plus one 'ignored' cluster!"); break;
				}
				t++;
			}
		}
		LMGroups.numGrpSeqs[n] = t;
		n++;
	}
	cout << endl;
    }
    LMGroups.numGroups = n;
    
    delete lmclusters;

} // end PhyloTree::readLikelihoodMappingGroups

//**************************************

void PhyloTree::doLikelihoodMapping() {
    // TODO For Heiko: Please add code here
    // vector<QuartetInfo> lmap_quartet_info;
    // vector<SeqQuartetInfo> lmap_seq_quartet_info;
    // int areacount[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // int cornercount[4] = {0, 0, 0, 0};
    int64_t resolved, partly, unresolved;
    int64_t qid;
    ofstream out;
    string filename;
    
    if(params->lmap_cluster_file != NULL) {
	// cout << "YYY: test reading" << params->lmap_cluster_file << endl;
        readLikelihoodMappingGroups(params->lmap_cluster_file, LMGroups);
    } else {
        LMGroups.numGroups = 0; /* no clusterfile -> un-initialized */
        int64_t recommended_quartets;
        if (aln->getNSeq() > 10) {
            recommended_quartets = (int64_t)25*aln->getNSeq();
        } else {
            int64_t n = aln->getNSeq();
            recommended_quartets = n*(n-1)*(n-2)*(n-3)/24;
        }
        if (params->lmap_num_quartets > 0 && params->lmap_num_quartets < recommended_quartets) {
            outWarning("Number of quartets is recommended to be at least " + convertInt64ToString(recommended_quartets) + " s.t. each sequence is sampled sufficiently");
        }

    }

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
    initsvg(svgout, LMGroups);

    string lmap_epsfilename = (string)params->out_prefix + ".lmap.eps";
    FILE *epsout;
    epsout = fopen(lmap_epsfilename.c_str(), "w");
    initeps(epsout, LMGroups);

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
        out << "SeqIDs\tlh1\tlh2\tlh3\tweight1\tweight2\tweight3\tarea\tcorner" << endl;
        // HAS - adding area7/corner3 output
        for (qid = 0; qid < params->lmap_num_quartets; qid++) {
            out << "(" << lmap_quartet_info[qid].seqID[0]+1 << ","
                << lmap_quartet_info[qid].seqID[1]+1 << ","
                << lmap_quartet_info[qid].seqID[2]+1 << ","
                << lmap_quartet_info[qid].seqID[3]+1 << ")"
                << "\t" << lmap_quartet_info[qid].logl[0] 
                << "\t" << lmap_quartet_info[qid].logl[1] 
                << "\t" << lmap_quartet_info[qid].logl[2]
                << "\t" << lmap_quartet_info[qid].qweight[0] 
                << "\t" << lmap_quartet_info[qid].qweight[1] 
                << "\t" << lmap_quartet_info[qid].qweight[2] // HAS - adding area7/corner3 output
                << "\t" << lmap_quartet_info[qid].area + 1
                << "\t" << lmap_quartet_info[qid].corner + 1 << endl;
        }

//        PhyloTree::reportLikelihoodMapping(out);

    }

    resolved   = areacount[0] + areacount[1] + areacount[2];
    partly     = areacount[3] + areacount[4] + areacount[5];
    unresolved = areacount[6];
	
#if 0  // for debugging purposes only (HAS)
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

#endif

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

} // end PhyloTree::doLikelihoodMapping




void PhyloTree::reportLikelihoodMapping(ofstream &out) {
    int64_t resolved, partly, unresolved;
    int64_t qid;
    size_t leafNum = PhyloTree::aln->getNSeq();
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

	out << "LIKELIHOOD MAPPING ANALYSIS" << endl;
	out << "---------------------------" << endl << endl;
	out << "Number of quartets: " << params->lmap_num_quartets;
	if (params->lmap_num_quartets < LMGroups.uniqueQuarts) {
		out << " (randomly chosen with replacement from "
			<< LMGroups.uniqueQuarts << " existing unique quartets)" << endl << endl;
	} 
	else {
		out << " (all unique quartets)" << endl << endl;
	}        
	out << "Quartet trees are based on the selected model of substitution." << endl << endl;

	if(LMGroups.numGroups == 1) {
		int n=0, t;
		out << "Sequences are not grouped in clusters. Using sequences:" << endl;
		for (t=0; t<LMGroups.numGrpSeqs[n]; t++){
		   out << "  " << LMGroups.GroupA[t]+1 << ". " 
			<< PhyloTree::aln->getSeqName(LMGroups.GroupA[t]) << endl;
		}
		out << endl << "Ordered as in user-given cluster file, numbers according to alignment order." << endl;
		out << "All other sequences have been ignored." << endl << endl;
	}
	if((LMGroups.numGroups > 1)&&(LMGroups.numGroups <= 4)) {
		int n, t;
		out << "Sequences are grouped into " << LMGroups.numGroups << " clusters." << endl << endl;
		for (n=0; n<LMGroups.numGroups; n++){
			out << "Cluster " << n+1 << " \"" << LMGroups.Name[n] << "\" lists " << LMGroups.numGrpSeqs[n] << " sequences: " << endl;
			for (t=0; t<LMGroups.numGrpSeqs[n]; t++){
			   switch(n){
				case 0: 
				   out << "  " << LMGroups.GroupA[t]+1 << ". " 
					<< PhyloTree::aln->getSeqName(LMGroups.GroupA[t]) << endl;
				   break;
				case 1: 
				   out << "  " << LMGroups.GroupB[t]+1 << ". " 
					<< PhyloTree::aln->getSeqName(LMGroups.GroupB[t]) << endl;
				   break;
				case 2: 
				   out << "  " << LMGroups.GroupC[t]+1 << ". " 
					<< PhyloTree::aln->getSeqName(LMGroups.GroupC[t]) << endl;
				   break;
				case 3: 
				   out << "  " << LMGroups.GroupD[t]+1 << ". " 
					<< PhyloTree::aln->getSeqName(LMGroups.GroupD[t]) << endl;
				   break;
				default: outError("Number of Likelihood Mapping groups too high! PLEASE report this to the developers!"); break;
	    		   }
			}
			out << endl;
		}
		out << "Ordered as in user-given cluster file, numbers according to alignment order." << endl;
		out << "All other sequences have been ignored." << endl << endl;
	}

        out << endl << endl;
	out << "LIKELIHOOD MAPPING STATISTICS" << endl;
	out << "-----------------------------" << endl << endl;

	switch(LMGroups.numGroups){
	   case 1: 
	      out << "           (a,b)-(c,d)                              (a,b)-(c,d)      " << endl; break;
	   case 2: 
	      out << "           (a,a)-(b,b)                              (a,a)-(b,b)      " << endl; break;
	   case 3: 
	      out << "           (a,b)-(c,c)                              (a,b)-(c,c)      " << endl; break;
	   case 4: 
	      out << "           (a,b)-(c,d)                              (a,b)-(c,d)      " << endl; break;
	   default: outError("Number of Likelihood Mapping groups too high! PLEASE report this to the developers!"); break;
	}
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
	switch(LMGroups.numGroups){
	   case 1: 
	      out << "(a,d)-(b,c)            (a,c)-(b,d)      (a,d)-(b,c)            (a,c)-(b,d) " 
	          << endl << endl; break;
	   case 2: 
	      out << "(a,b)-(a,b)            (a,b)-(a,b)      (a,b)-(a,b)            (a,b)-(a,b) " 
	          << endl << endl; break;
	   case 3: 
	      out << "(a,c)-(b,c)            (a,c)-(b,c)      (a,c)-(b,c)            (a,c)-(b,c) " 
	          << endl << endl; break;
	   case 4: 
	      out << "(a,d)-(b,c)            (a,c)-(b,d)      (a,d)-(b,c)            (a,c)-(b,d) " 
	          << endl << endl; break;
	   default: outError("Number of Likelihood Mapping groups too high! PLEASE report this to the developers!"); break;
	}
        //      |<---                                   80 chars                           --->|
	out << "Division of the likelihood mapping plots into 3 or 7 areas." << endl;
	out << "On the left the areas show support for one of the different groupings" << endl;
	out << "like (a,b|c,d)." << endl;

	out << "On the right the right quartets falling into the areas 1, 2 and 3 are" << endl;
	out << "informative. Those in the rectangles 4, 5 and 6 are partly informative" << endl;
	out << "and those in the center (7) are not informative." << endl << endl;

	switch(LMGroups.numGroups){
	   case 1: 
		out << "Sequences a,b,c,d are drawn from all included sequences." << endl << endl; break;
	   case 2: 
		out << "Sequences a(2x) and b(2x) are drawn from clusters 1 and 2, respectively, with" << endl;
		out << "Cluster 1: " << LMGroups.Name[0] << endl;
		out << "Cluster 2: " << LMGroups.Name[1] << endl << endl;
		break;
	   case 3: 
		out << "Sequences a, b and c(2x) are drawn from clusters 1, 2 and 3, respectively, with" << endl;
		out << "Cluster 1: " << LMGroups.Name[0] << endl;
		out << "Cluster 2: " << LMGroups.Name[1] << endl;
		out << "Cluster 3: " << LMGroups.Name[2] << endl << endl;
		break;
	   case 4: 
		out << "Sequences a,b,c,d are drawn from clusters 1, 2, 3 and 4, respectively, with" << endl;
		out << "Cluster 1: " << LMGroups.Name[0] << endl;
		out << "Cluster 2: " << LMGroups.Name[1] << endl;
		out << "Cluster 3: " << LMGroups.Name[2] << endl;
		out << "Cluster 4: " << LMGroups.Name[3] << endl << endl;
		break;
	   default: outError("Number of Likelihood Mapping groups too high! PLEASE report this to the developers!"); break;
	}

	out << "Note, that the corners only make a difference if the sequences are" << endl;
	out << "clustered in groups. Furthermore, while sequences should occur about" << endl;
	out << "equally often in unclustered mappings, in clustered mappings their" << endl;
	out << "occurrence rates depend on the group sizes the quartets are drawn from." << endl << endl;

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
	    int64_t sumq0, sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];
	    if (sumq>0) sumq0=sumq;
	    else sumq0=1;

	    if (qid < leafNum) {
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << setw(4) << qid+1 
                   << setw(9) << sumq
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[7]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[7]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[8]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[8]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[9]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[9]/sumq0 << ")  " 
        	   << PhyloTree::aln->getSeqName(qid) << endl;
	    } else {
	       out << "-----------------------------------------------------------------------------" << endl;
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << "    "
                   << setw(9) << sumq
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[7]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[7]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[8]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[8]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[9]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[9]/sumq0 << ")  " << endl;
	    }
        }

        out << endl << endl << "Quartet support of areas 1-7 (mainly for clustered analysis):" << endl << endl;
        out << "                   resolved                                           partly                                             unresolved  name" << endl;
        out << "     #quartets     1 (% 1)          2 (% 2)          3 (% 3)          4 (% 4)          5 (% 5)          6 (% 6)          7 (% 7)" << endl;
	out << "------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
        for (qid = 0; qid <= leafNum; qid++) {
	    int64_t sumq0, sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];
	    if (sumq>0) sumq0=sumq;
	    else sumq0=1;

	    if (qid < leafNum) {
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << setw(4) << qid+1 
                   << setw(9) << sumq
                   << setw(7) << lmap_seq_quartet_info[qid].countarr[0] 
		   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[0]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[1] 
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[1]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[2]
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[2]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[3]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[3]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[4]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[4]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[5]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[5]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[6]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[6]/sumq0 << ")  "
        	   << PhyloTree::aln->getSeqName(qid) << endl;
	    } else {
	       out << "------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << "    "
                   << setw(9) << sumq
                   << setw(7) << lmap_seq_quartet_info[qid].countarr[0] 
		   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[0]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[1] 
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[1]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[2]
                   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[2]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[3]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[3]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[4]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[4]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[5]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[5]/sumq0 << ") "
        	   << setw(7) << lmap_seq_quartet_info[qid].countarr[6]
        	   << " (" << setw(6) << (double) 100.0*lmap_seq_quartet_info[qid].countarr[6]/sumq0 << ")  " 
	           << endl << endl;
	    }
        }

        out << endl << endl << "Quartet resolution per sequence (phylogenetic information):" << endl << endl;
        out << "      #quartets   resolved           partly          unresolved  name" << endl;
	out << "-----------------------------------------------------------------------------" << endl;
        for (qid = 0; qid <= leafNum; qid++) {
	    int64_t resolved = lmap_seq_quartet_info[qid].countarr[LM_REG1] + lmap_seq_quartet_info[qid].countarr[LM_REG2] + lmap_seq_quartet_info[qid].countarr[LM_REG3];
	    int64_t partly   = lmap_seq_quartet_info[qid].countarr[LM_REG4] + lmap_seq_quartet_info[qid].countarr[LM_REG5] + lmap_seq_quartet_info[qid].countarr[LM_REG6];
	    int64_t unres    = lmap_seq_quartet_info[qid].countarr[LM_REG7];
	    int64_t sumq0, sumq = lmap_seq_quartet_info[qid].countarr[7] + lmap_seq_quartet_info[qid].countarr[8] + lmap_seq_quartet_info[qid].countarr[9];
	    if (sumq>0) sumq0=sumq;
	    else sumq0=1;

	    if (qid < leafNum) {
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << setw(4) << qid+1 
                   << setw(9) << sumq
                   << setw(7) << resolved
		   << " (" << setw(6) << (double) 100.0*resolved/sumq0 << ") "
        	   << setw(7) << partly
                   << " (" << setw(6) << (double) 100.0*partly/sumq0 << ") "
        	   << setw(7) << unres
                   << " (" << setw(6) << (double) 100.0*unres/sumq0 << ")  "
        	   << PhyloTree::aln->getSeqName(qid) << endl;
	    } else {
	       out << "-----------------------------------------------------------------------------" << endl;
               out.setf(ios::fixed, ios::floatfield); // set fixed floating format
               out.precision(2);
               out << "    "
                   << setw(9) << sumq
                   << setw(7) << resolved
		   << " (" << setw(6) << (double) 100.0*resolved/sumq0 << ") "
        	   << setw(7) << partly
                   << " (" << setw(6) << (double) 100.0*partly/sumq0 << ") "
        	   << setw(7) << unres
                   << " (" << setw(6) << (double) 100.0*unres/sumq0 << ")  " << endl;
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

