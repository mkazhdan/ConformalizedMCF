<html><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>Conformalized MCF</title>
<STYLE>
DL{ margin: 0px 0; }
UL{ list-style: none;}
</STYLE>
</head>
<body onload="_init();">
<center><h1>Can Mean-Curvature Flow be Modified to be Non-singular?<br> (Version 2.0)</h1></center>
<center>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version2/#LINKS">links</a>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version2/#EXECUTABLES">executables</a>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version2/#COMPILING">compiling</a>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version2/#CHANGES">changes</a>
</center>
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<a href="http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf">Kazhdan, 2012 Paper</a>
<a href="https://dl.acm.org/citation.cfm?id=3073615">Aigerman et al., 2017</a><br>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version2/ConformalizedMCF.exe.zip">Win64 Executables</a><br>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version2/ConformalizedMCF.zip">Source Code</a><br>
<a href="https://github.com/mkazhdan/ConformalizedMCF">GitHub Repository</a><br>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version2/license.txt">License</a><br>
<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>

<UL>
	<li><font size="+1"><b>ConformalizedMCF</b></font>
		<UL>
			<li>
				<DL>
					<DT><U>Command line arguments</U>
					<DD>
					<DL>
						<dt><b>--in</b> &lt;<i>input mesh</i>&gt;
						<dd> This string is the name of the file containing the triangulated model.
						It is assumed that the file is in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and the file-name ends with a ".ply" extension.

						<dt>[<b>--flow</b> &lt;<i>type of surface flow</i>]&gt;
						<dd> This argument specifies the type of flow used to evolve the surface. "<b>1</b>" corresponds to traditional mean-curvature flow. "<b>2</b>" corresponds to the conformalized mean-curvature flow. And "<b>3</b>" corresponds to heat flow.<br>
						If no value specified, a default value of "<b>2</b>" is used.

						<dt>[<b>--outHeader</b> &lt;<i>output file headers</i>&gt;]
						<dd> This string specifies the header for the filenames used when outputting the evolving surface throughout the course of the flow.<br>
						If no argument is specified, no output files are generated.

						<dt>[<b>--steps</b> &lt;<i>number of steps of evolution to be performed</i>&gt;]
						<dd> This argument specifies the number of anisotropy normalizations that are to be preformed. The parameter for this argument is specified as either <i>start:increment:end</i>, <i>start:end</i>, or <i>start</i> where <i>start</i> is the first evolution step at which files are output, <i>end</i> is the number of evolution steps performed, and <i>increment</i> is the number of evolution steps between file outputs. If <i>increment</i> is not set, a default value of 1 is used. If <i>end</i> is not set, the value is set to <i>start</i>.<br>
						If no argument is specified, no evolution is performed.

						<dt>[<b>--stepSize</b> &lt;<i>temporal discretization</i>&gt;]
						<dd> This argument specifies the size of the temporal discretization (with smaller values resulting in slower flows and larger values resulting in faster ones).<br>
						If not argument is specified, a default value of 0.0001 is used.

						<dt> [<b>--verbose</b>]
						<dd> It this flag is specified error measures are computed at each step of the evolution and output to STDOUT.
					</DL>
			</DL>
		</UL>

   <li><font size="+1"><b>ConformalizedMCFOrbifoldVisualization</b></font>
		<UL>
			<li>
				<DL>
					<DT><U>Command line arguments</U>
					<DD>
					<dl>
						<dt><b>--in</b> &lt;<i>input mesh</i>&gt;
						</dt><dd> This string is the name of the file containing the triangulated model.
						It is assumed that the file is in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and the file-name ends with a ".ply" extension.<br>
						The model must have either disk or sphere topology.

						<dt>[<b>--sym</b> &lt;<i>type of orbifold symmetry</i>]&gt;
						<dd> This argument specifies the type of orbifold symmetry to be used. Options are:
						<ul>
							<li> <b>C</b>&lt;<i>N</i>&gt;: Cyclic symmetry group of order <i>N</i>
							<li> <b>D</b>&lt;<i>N</i>&gt;: Dihedral symmetry group of order <i>2N</i>
							<li> <b>T</b>: Tetrahedral symmetry group
							<li> <b>O</b>: Octahedral symmetry group
							<li> <b>I</b>: Icosahedral symmetry group
						</ul>
						If no value specified, a default value of "<b>C1</b>" is used.

						<dt>[<b>--stepSize</b> &lt;<i>temporal discretization</i>&gt;]
						<dd> This argument specifies the size of the temporal discretization (with smaller values resulting in slower flows and larger values resulting in faster ones).<br>
						If not argument is specified, a default value of 1.0 is used.

						<dt> [<b>--verbose</b>]
						<dd> It this flag is specified error measures are computed at each step of the evolution and output to STDOUT.
					</dl>					
			</DL>
				<DL>
					<DT><U>Viewer controls</U>
					<DD>
						'<I><B>q</B></I>', '<I><B>w</B></I>', '<I><B>a</B></I>', '<I><B>s</B></I>', '<I><B>z</B></I>', '<I><B>x</B></I>': rotate<BR>
						<B>[UP ARROW]</B>, <B>[DOWN ARROW]</B>: zoom<BR>
						<B>[SPACE]</B>: Toggle flow<BR>
						See bottom right panel for additional key bindings
			</DL>
		</UL>
</ul>


<hr>
<a name="COMPILING"><b>SOURCE CODE COMPILATION</b></a><br>
The source code requires the use of a numerical solver. By default, the <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a> solver is used. As the default implementation may be slow, the code supports other solvers. Specifically, if you have the <a href="https://software.intel.com/en-us/mkl">Intel Math Kernel Library</a> installed, you can define <b>EIGEN_USE_MKL_ALL</b> to use a faster solver. Alternatively, if you have the <a href="http://www.cise.ufl.edu/research/sparse/cholmod/">CHOLMOD</a> library installed, you can use that by defineing <b>USE_CHOLMOD</b> (and also <b>USE_SUITESPARSE</b> if you are using the <a href="http://faculty.cse.tamu.edu/davis/suitesparse.html">SuiteSparse</a> version).

<hr>
<a name="CHANGES"><b>CHANGES</b></a><br>
<a href="http://www.cs.jhu.edu/~misha/Code/ConformalizedMCF/Version1/">Version 1.0</a>:
<ol>
<li> Added spherical orbifold visualization code.
</li></ol>


<hr>
<a href="http://www.cs.jhu.edu/~misha">HOME</a>


</body></html>