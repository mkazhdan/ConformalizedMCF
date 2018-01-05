
<html>
<head>
<title>Conformalized MCF</title>
</head>
<body ONLOAD="_init();">
<CENTER><H1>Can Mean-Curvature Flow be Modified to be Non-singular?<BR> (Version 2.0)</A></H1></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#EXECUTABLE">executable</A>
<A HREF="#COMPILING">compiling</A>
<A HREF="#CHANGES">changes</A>
</CENTER>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf">Kazhdan, 2012 Paper</A>
<A HREF="https://dl.acm.org/citation.cfm?id=3073615">Aigerman et al., 2017</A><BR>
<A HREF="ConformalizedMCF.exe.zip">Win64 Executables</A></br>
<A href="ConformalizedMCF.zip">Source Code</A><br>
<A href="license.txt">License</A><br>
<HR>
<A NAME="EXECUTABLE"><B>EXECUTABLE</B></A><br>

<UL>
<DL>
<FONT SIZE="+1" ><B>ConformalizedMCF</B></FONT>
<DT><b>--in</b> &#60;<i>input mesh</i>&#62;
<DD> This string is the name of the file containing the triangulated model.
It is assumed that the file is in the
<A HREF="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format and the file-name ends with a ".ply" extension.

<DT>[<b>--flow</b> &#60;<i>type of surface flow</i>]&#62;
<DD> This argument specifies the type of flow used to evolve the surface. "<B>1</B>" corresponds to traditional mean-curvature flow. "<B>2</B>" corresponds to the conformalized mean-curvature flow. And "<B>3</B>" corresponds to heat flow.<BR>
If no value specified, a default value of "<B>2</B>" is used.

<DT>[<b>--outHeader</B> &#60;<i>output file headers</i>&#62;]
<DD> This string specifies the header for the filenames used when outputting the evolving surface throughout the course of the flow.<BR>
If no argument is specified, no output files are generated.

<DT>[<b>--steps</b> &#60;<i>number of steps of evolution to be performed</i>&#62;]
<DD> This argument specifies the number of anisotropy normalizations that are to be preformed. The parameter for this argument is specified as either <I>start:increment:end</I>, <I>start:end</I>, or <I>start</I> where <I>start</I> is the first evolution step at which files are output, <I>end</I> is the number of evolution steps performed, and <I>increment</I> is the number of evolution steps between file outputs. If <I>increment</I> is not set, a default value of 1 is used. If <I>end</I> is not set, the value is set to <I>start</I>.<BR>
If no argument is specified, no evolution is performed.

<DT>[<b>--stepSize</b> &#60;<i>temporal discretization</i>&#62;]
<DD> This argument specifies the size of the temporal discretization (with smaller values resulting in slower flows and larger values resulting in faster ones).<BR>
If not argument is specified, a default value of 0.0001 is used.

<DT> [<B>--verbose</B>]
<DD> It this flag is specified error measures are computed at each step of the evolution and output to STDOUT.

</DL>
</UL>

<UL>
<DL>
<FONT SIZE="+1" ><B>ConformalizedMCFOrbifoldVisualization</B></FONT>

<DT><b>--in</b> &#60;<i>input mesh</i>&#62;
<DD> This string is the name of the file containing the triangulated model.
It is assumed that the file is in the <A HREF="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format and the file-name ends with a ".ply" extension.<br>
The model must have either disk or sphere topology.

<DT>[<b>--sym</b> &#60;<i>type of orbifold symmetry</i>]&#62;
<DD> This argument specifies the type of orbifold symmetry to be used. Options are:
<UL>
<LI> <B>C</B>&#60;<I>N</I>&#62;: Cyclic symmetry group of order <I>N</I>
<LI> <B>D</B>&#60;<I>N</I>&#62;: Dihedral symmetry group of order <I>2N</I>
<LI> <B>T</B>: Tetrahedral symmetry group
<LI> <B>O</B>: Octahedral symmetry group
<LI> <B>I</B>: Icosahedral symmetry group
</UL>
If no value specified, a default value of "<B>C1</B>" is used.

<DT>[<b>--stepSize</b> &#60;<i>temporal discretization</i>&#62;]
<DD> This argument specifies the size of the temporal discretization (with smaller values resulting in slower flows and larger values resulting in faster ones).<BR>
If not argument is specified, a default value of 1.0 is used.

<DT> [<B>--verbose</B>]
<DD> It this flag is specified error measures are computed at each step of the evolution and output to STDOUT.

</DL>
</UL>


<HR>
<A NAME="COMPILING"><B>SOURCE CODE COMPILATION</B></A><br>
The source code requires the use of a numerical solver. By default, the <A HREF="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</A> solver is used. As the default implementation may be slow, the code supports other solvers. Specifically, if you have the <A HREF="https://software.intel.com/en-us/mkl">Intel Math Kernel Library</A> installed, you can define <B>EIGEN_USE_MKL_ALL</B> to use a faster solver. Alternatively, if you have the <A HREF="http://www.cise.ufl.edu/research/sparse/cholmod/">CHOLMOD</A> library installed, you can use that by defineing <B>USE_CHOLMOD</B> (and also <B>USE_SUITESPARSE</B> if you are using the <A HREF="http://faculty.cse.tamu.edu/davis/suitesparse.html">SuiteSparse</A> version).

<HR>
<A NAME="CHANGES"><B>CHANGES</B></A><br>
<A HREF="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version1/">Version 1.0</A>:
<OL>
<LI> Added spherical orbifold visualization code.
</OL>


<HR>
<A HREF="http://www.cs.jhu.edu/~misha">HOME</A>
</body>
</html>
