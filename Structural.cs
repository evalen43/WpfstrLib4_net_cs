using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data.SQLite;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using System.Xml;
using System.Collections;
using MatrixCs;

namespace WpfstrLib4_net_cs
{
    public class Structural
    {
        public struct Node
        {
            public string nodeid;
            public double x, y, z;
            public double distance(Node n)
            {
                double length = Math.Sqrt(Math.Pow((this.x - n.x), 2) + Math.Pow((this.y - n.y), 2) + Math.Pow((this.z - n.z), 2));
                return length;
            }
        }
        public struct Member
        {
            public string memberid;
            public int sec_id, mat_id, inc1, inc2;
            public double d, cx, cy, cz, beta;
        }
        public struct Section
        {
            public string sectionid, sectiontype;
            public double W, A, Ix, Zx, Sx, rx, Iy, Zy, Sy, ry, J, beta, od, wth;
            public void SecTube()
            {
                double id = od - 2 * wth;
                double area = 0.25 * Math.PI * (od * od - id * id);
                double iner = Math.PI * (od * od * od * od - id * id * id * id) / 64;
                double s = 2 * iner / od;
                double r = Math.Sqrt(iner / area);
                A = area;
                A = area;
                Ix = iner;
                Iy = iner;
                J = 2 * iner;
                Sx = s;
                Zx = (Math.Pow(od, 3) - Math.Pow(id, 3)) / 6;
                Zy = Zx;
                Sy = s;
                rx = r;
                ry = r;
            }
            public void Rectangular()
            {
                double b = wth;
                double d = od;
                double area = d * b;
                double inerx = b * Math.Pow(d, 3) / 12;
                double inery = d * Math.Pow(b, 3) / 12;
                rx = Math.Sqrt(inerx / area);
                ry = Math.Sqrt(inery / area);
                double beta = (1 / 3) - 0.21 * (b / d) * (1 - Math.Pow(b / d, 4) / 12);
                J = beta * Math.Pow(b, 3) * d;
                A = area;
                Ix = inerx;
                Iy = inery;//  iner;
                Sx = Ix * 2 / d;
                Sy = Ix * 2 / b;

            }

            public bool AISC(double scaleL, double scaleF)
            {
                string cs = @"URI=file:C:\Git\repos\WpfstrLib4_net_cs\aisc_shapes_v15_US.db3";
                using var db = new SQLiteConnection(cs);
                db.Open();
                string stm = "SELECT * FROM shapesv15_US"; /// name of table in database (shapesv15_US)
                SQLiteCommand cmd = new(stm, db);
                SQLiteDataReader rdr = cmd.ExecuteReader();

                bool found = false;
                while (rdr.Read()) /// Loop to execute the query step by step, to get rows of result
                {
                    string? mystring = rdr[0].ToString();
                    if (mystring == sectionid)
                    {
                        found = true;
                        W = Convert.ToDouble(rdr[1].ToString()) * scaleL / scaleF;       ///linear weight in lbf/ft before conversion to kN/m
                        A = Convert.ToDouble(rdr[2].ToString()) * Math.Pow(scaleL, 2);  /// area in inch4 before conversion to m2
                        Ix = Convert.ToDouble(rdr[3].ToString()) * Math.Pow(scaleL, 4); /// strong axis moment of inertia in inch4 before conversion to m4
                        Zx = Convert.ToDouble(rdr[4].ToString()) * Math.Pow(scaleL, 3); ///strong axis plastic module in inch3 before conversion to m3
                        Sx = Convert.ToDouble(rdr[5].ToString()) * Math.Pow(scaleL, 3); /// strong axis elastic module in inch3 before conversion to m3
                        rx = Convert.ToDouble(rdr[6].ToString()) * scaleL;              /// strong axis radius of gyration in inches before conversion to m
                        Iy = Convert.ToDouble(rdr[7].ToString()) * Math.Pow(scaleL, 4); /// weak axis moment of inertia in inch4 before conversion to m4
                        Zy = Convert.ToDouble(rdr[8].ToString()) * Math.Pow(scaleL, 3); /// weak axis plastic module in inch3 before conversion to m3
                        Sy = Convert.ToDouble(rdr[9].ToString()) * Math.Pow(scaleL, 3); /// weak axis elastic module in inch3 before conversion to m3
                        ry = Convert.ToDouble(rdr[10].ToString()) * scaleL;              /// weak axis radius of gyration in inch before conversion to m
                        J = Convert.ToDouble(rdr[11].ToString()) * Math.Pow(scaleL, 4); /// polar moment of inertia in inch4 before conversion to m4
                        break;
                    }
                }
                return found;
            }

        }
        public struct Material
        {
            public int matid;
            public string matname, mtablename;
            public double Ematerial, poisson, denMaterial, Gmod;
        }
        public struct Boundary
        {
            public int[] ib;//= new int[4] { -1, 0, 0, 0 };
            public string nodebid;// = "N/A";
            public string nbtype;// = "ENCASTRE";
            public double[] w;

            public void SetBound(string strutype)
            {

                //ib![0] = int.Parse(nodebid);
                if (strutype == "Frame2D")
                {
                    if (nbtype == "ENCASTRE")
                    {
                        ib[1] = 0; ib[2] = 0; ib[3] = 0;
                        w[0] = 0; w[1] = 0.0; w[2] = 0.0;
                    }
                    else if (nbtype == " HINGE")
                    {
                        ib[1] = 0; ib[2] = 0; ib[3] = 1;
                        w[0] = 0; w[1] = 0.0;
                    }
                    else if (nbtype == "EFSLIDEX")
                    {
                        ib[1] = 1; ib[2] = 0; ib[3] = 0;
                    }
                    else if (nbtype == "EFSLIDEY")
                    {
                        ib[1] = 0; ib[2] = 1; ib[3] = 0;
                    }
                }

            }

        }
        public struct Jload
        {
            public double[] P;
            public int jnodeid, jlc;// = 0, jlc = 0;
            public string nodeid;// = "";
        }
        public struct Mload
        {
            public double[] f;
            public double[] vlocal;
            public string mem_id, loadtype, sysref;
            public int lc_number, elemindex;// = 0, elemindex = 0;
            public double wa, wb, a, P;// = 0, wb = 0, a = 0, P = 0;
        }

        public class StructuralCs
        {


            public Node[] wxcoor;
            public Member[] wxelement;
            public Material[] wxmaterial;
            public Section[] wxsection;
            public Boundary[] wxboundary;
            public Jload[] wjload;
            public Mload[] wmload;
            public string[] nlist, nelemlist, nmatlist, nseclist, nbnlist, njload, nmldlist;

            public ArrayList wxjload = new();
            public ArrayList wxmload = new();

            private string unitL, unitF, unitS, structure, code, codetext, project, unitDen, caseno, loadtype;
            private int ndf, nne, nn, nbn, ne, nlc, ndfel, nmat, n, nsec, ms;
            private double[] intforc;

            private double meters, scaleL, scaleF, scaleS, scaleDen, fyield;
            private readonly string filein, fileout, filelog;
            private readonly double g = 9.806;
            private double[,] tk, al, reac, fem_dload;
            //private Matrix[,] tk, al, reac,fem_dload;
            public string Structure { get => structure; set => structure = value; }
            public string Project { get => project; set => project = value; }
            public string Codetext { get => codetext; set => codetext = value; }
            public int Ndfel { get => ndfel; set => ndfel = value; }
            public int Ne { get => ne; set => ne = value; }
            public int Nn { get => nn; set => nn = value; }
            public int N { get => n; set => n = value; }
            public int Ms { get => ms; set => ms = value; }

            public StructuralCs(string filein, string fileout, string filelog)
            { this.filein = filein; this.fileout = fileout; this.filelog = filelog; }
            public void Assem()
            {
                //tk = new double[n, ms];
                for (int i = 0; i < ne; i++)
                {
                    //Matrix rot = Rotmat(i);
                    //Matrix elst = elem_stiff(i);
                    double[,] rot = Rotmat(i);
                    double[,] elst = elem_stiff(i);
                    if (i == 0)
                    {
                        MatPrint(rot, "Rotation Matrix", false);
                        MatPrint(elst, "local stiffness elst", true);
                    }
                    //elst = ~rot * elst * rot;//~rot is rot transposed
                    elst = BTAB3(elst, rot);
                    if (i == 0) MatPrint(elst, "global stiffness elst", true);
                    elassgen(i, elst);
                }
#if DEBUG
                MatPrint(tk, "Stiffness Matrix after Assembly", true);
                //using (StreamWriter Debugfile = new(filelog, append: false))
                //{
                //    Debugfile.WriteLine("Stiffness Matrix after Assembly");
                //    for (int i = 0; i < n; i++)
                //    {
                //        for (int j = 0; j < ms; j++)
                //        {
                //            Debugfile.Write(tk![i, j].ToString() + "\t");
                //        }
                //        Debugfile.WriteLine("");
                //    }
                //}
#endif
            }
            public double[,] Rotmat(int nel)
            {
                /// <summary>
                /// The rotation matrix represents the relationship between local (x',y') and global (x,y) coordinates
                /// x'=cos(ang)*x+sin(ang)*y
                /// y'=-sin(ang)*x+cos(ang)*y
                /// </summary>
                /// <param name=nel></param>
                /// <returns>rotation matrix</returns>
                double CXZ, COSBET, SINBET, beta;
                double cx, cy, cz;
                //Matrix rot = new(ndfel, ndfel);

                double[,] rot = new double[ndfel, ndfel];
                cx = wxelement[nel].cx;// c[0];
                cy = wxelement[nel].cy;// c[1];
                cz = wxelement[nel].cz;// c[2];
                beta = wxelement[nel].beta;// c[3];
                CXZ = Math.Sqrt(cx * cx + cz * cz);
                switch (structure)
                {
                    case "Truss2D":
                        {
                            // ----
                            // 2D TRUSS
                            // ----
                            rot[0, 0] = cx;
                            rot[0, 1] = cy;
                            rot[1, 0] = -cy;
                            rot[1, 1] = cx;
                            break;
                        }

                    case "Truss3D":
                        {
                            // ----
                            // Space TRUSS
                            // ----
                            if ((CXZ != 0.0))
                            {
                                rot[0, 0] = cx;
                                rot[0, 1] = cy;
                                rot[0, 2] = cz;
                                rot[1, 0] = -cx * cy / CXZ;
                                rot[1, 1] = CXZ;
                                rot[1, 2] = -cy * cz / CXZ;
                                rot[2, 0] = -cz / CXZ;
                                rot[2, 1] = 0.0;
                                rot[2, 2] = cx / CXZ;
                            }
                            else if ((CXZ == 0))
                            {
                                rot[0, 0] = 0.0;
                                rot[0, 1] = cy;
                                rot[0, 2] = 0.0;
                                rot[1, 0] = -cy;
                                rot[1, 1] = 0.0;
                                rot[1, 2] = 0.0;
                                rot[2, 0] = 0.0;
                                rot[2, 1] = 0.0;
                                rot[2, 2] = 1.0;
                            }
                            break;
                        }
                    case "Frame2D":
                        {
                            // ------------
                            // PLANE FRAMES
                            // ------------
                            rot[0, 0] = cx;
                            rot[0, 1] = cy;
                            rot[1, 0] = -cy;
                            rot[1, 1] = cx;
                            rot[2, 2] = 1.0;
                            break;
                        }
                    case "Grid":
                        {
                            // ------------
                            // Grids
                            // ------------
                            rot[0, 0] = cx;
                            rot[0, 1] = cy;
                            rot[1, 0] = -cy;
                            rot[1, 1] = cx;
                            rot[2, 2] = 1.0;
                            break;
                        }
                    case "Frame3D":
                        {
                            // ----
                            // Space Frames
                            // ----
                            COSBET = Math.Cos(beta);
                            SINBET = Math.Sin(beta);
                            if ((CXZ != 0.0))
                            {
                                rot[0, 0] = cx;
                                rot[0, 1] = cy;
                                rot[0, 2] = cz;
                                rot[1, 0] = (-cx * cy * COSBET - cz * SINBET) / CXZ;
                                rot[1, 1] = CXZ * COSBET;
                                rot[1, 2] = (-cy * cz * COSBET + cx * SINBET) / CXZ;
                                rot[2, 0] = (cx * cy * SINBET - cz * COSBET) / CXZ;
                                rot[2, 1] = -CXZ * SINBET;
                                rot[2, 2] = (cy * cz * SINBET + cx * COSBET) / CXZ;
                            }
                            else if ((CXZ == 0.0))
                            {
                                // ----------------------------------------------
                                // THE MEMBER Is VERTICAL And PERPENDICULAR TO XZ
                                // ----------------------------------------------
                                rot[0, 1] = 1.0;
                                rot[1, 0] = -COSBET;
                                rot[1, 2] = SINBET;
                                rot[2, 0] = SINBET;
                                rot[2, 2] = COSBET;
                            }
                            break;
                        }
                }
                for (int I = 0; I < ndf; I++)
                {
                    for (int J = 0; J < ndf; J++)
                        rot[I + ndf, J + ndf] = rot[I, J];
                }
                return rot;
            }
            public void MatPrint(double[,] A, string label, bool xappend)
            {
                int irow, icol;
                irow = A.GetUpperBound(0) + 1;
                icol = A.GetUpperBound(1) + 1;

                using (StreamWriter Debugfile = new(filelog, append: xappend))
                {

                    Debugfile.WriteLine(label);
                    Debugfile.WriteLine("Nrows " + irow.ToString());
                    Debugfile.WriteLine("Ncols " + icol.ToString());
                    for (int i = 0; i < irow; i++)
                    {
                        for (int j = 0; j < icol; j++)
                            Debugfile.Write(A[i, j].ToString() + "\t");
                        Debugfile.WriteLine("");
                    }
                    Debugfile.Close();
                }
            }
            public void MatrixPrint(Matrix A, string label, int ndfel, bool xappend)
            {
                int irow, icol;
                irow = ndfel;// A[]GetUpperBound(0);
                icol = ndfel;// A.GetUpperBound(1);

                using (StreamWriter Debugfile = new StreamWriter(filelog, append: xappend))
                {

                    Debugfile.WriteLine(label);
                    Debugfile.WriteLine("Nrows " + irow.ToString());
                    Debugfile.WriteLine("Ncols " + icol.ToString());
                    for (int i = 0; i < irow; i++)
                    {
                        for (int j = 0; j < icol; j++)
                            Debugfile.Write(A[i, j].ToString() + "\t");
                        Debugfile.WriteLine("");
                    }
                    Debugfile.Close();
                }
            }
            public void Bandgauss()
            {
                int n1, k1, l, ni, k2, k3, k;
                double c;
                double[] d = new double[n];
                n1 = n - 1;
                for (k = 1; k <= n1; k++)
                {
                    c = tk![k - 1, 0];
                    k1 = k + 1;
                    if (Math.Abs(c) <= 1.0e-06) //l_2:format(1ho, "**** singularity in row", i5, 1x, "****")
                        throw new DivideByZeroException();
                    ni = k1 + ms - 2;
                    l = Math.Min(ni, n);
                    for (int j = 2; j <= ms; j++)
                        d[j - 1] = tk[k - 1, j - 1];

                    for (int j = k1; j <= l; j++)
                    {
                        k2 = j - k + 1;
                        tk[k - 1, k2 - 1] = tk[k - 1, k2 - 1] / c;  // divide row by diagonal coefficient
                    }
                    for (int icol = 1; icol <= nlc; icol++)
                        al![k - 1, icol - 1] = al[k - 1, icol - 1] / c;
                    for (int i = k1; i <= l; i++)
                    {
                        k2 = i - k1 + 2;
                        c = d[k2 - 1];
                        for (int j = i; j <= l; j++)
                        {
                            k2 = j - i + 1;
                            k3 = j - k + 1;
                            tk[i - 1, k2 - 1] = tk[i - 1, k2 - 1] - c * tk[k - 1, k3 - 1];
                        }
                        for (int icol = 1; icol <= nlc; icol++)
                            al[i - 1, icol - 1] = al[i - 1, icol - 1] - c * al[k - 1, icol - 1]; // eliminate unknown x(k) from row i
                    }
                }
                ///compute last unknown
                if (Math.Abs(tk[n - 1, 0]) < 1.0e-06) //l_17: write(iout, 2) n
                    throw new DivideByZeroException();

                for (int icol = 1; icol <= nlc; icol++)
                    al[n - 1, icol - 1] = al[n - 1, icol - 1] / tk[n - 1, 1 - 1];

                ///apply backsubstitution process to compute remaining unknowns
                for (int i = 1; i <= n1; i++)
                {
                    k = n - i;
                    k1 = k + 1;
                    ni = k1 + ms - 2;
                    l = Math.Min(ni, n);
                    for (int j = k1; j <= l; j++)
                    {
                        k2 = j - k + 1;
                        for (int icol = 1; icol <= nlc; icol++)
                            al![k - 1, icol - 1] = al[k - 1, icol - 1] - tk[k - 1, k2 - 1] * al[j - 1, icol - 1];
                    }
                }
            } // end
            public void InputXML()
            {
                project = " ";
                codetext = " ";
                char[] delims = new[] { '\r', '\n' };
                XmlDocument xDoc = new();
                xDoc.Load(filein);
                XmlElement root = xDoc.DocumentElement;
                structure = root.Name;
                switch (structure)
                {
                    case ("Frame2D"):
                        ndf = 3; nne = 2; ndfel = ndf * nne;
                        break;
                    case ("Frame3D"):
                        ndf = 6; nne = 2; ndfel = ndf * nne;
                        break;
                    case ("Trust3D"):
                        ndf = 3; nne = 2; ndfel = ndf * nne;
                        break;
                    case ("Trust2D"):
                        ndf = 2; nne = 2; ndfel = ndf * nne;
                        break;
                    default:
                        ndf = 2; nne = 2; ndfel = ndf * nne;
                        break;
                    case ("Grid"):
                        ndf = 3; nne = 2; ndfel = ndf * nne;
                        break;
                }
                foreach (XmlNode node in xDoc.DocumentElement!.ChildNodes)
                {
                    Section nsection = new();
                    if (node.Name == "project") project = node.InnerText.ToString().Trim();
                    else if (node.Name == "code")
                    {
                        if (node.Attributes!.Count > 0) { unitS = node.Attributes[0].Value; }
                        else unitS = "kN/m2";
                        scaleS = TokNperM2(unitS);
                        codetext = node.InnerText.ToString().Trim();
                        string[] words = codetext.Split(new char[] { ' ', '=' });
                        for (int k = 0; k < words.Length; k++)
                        {
                            if (words[k] == "code") code = words[k + 1];
                            if (words[k] == "fy") fyield = Convert.ToDouble(words[k + 1]) * scaleS;
                        }
                        Debug.WriteLine("Code" + "\t" + code + "\t" + "fyield" + "\t" + fyield.ToString());
                    }
                    else if (node.Name == "nodes")
                    {
                        if (node.Attributes!.Count > 0) { unitL = node.Attributes[0].Value; }
                        else unitL = "m";
                        scaleL = ToMeters(unitL);
                        string content = node.InnerText.ToString().Trim();
                        string[] lines = content.Split(delims, StringSplitOptions.RemoveEmptyEntries);
                        nn = lines.Length;
                        n = nn * ndf;
                        wxcoor = new Node[nn];
                        nlist = new string[nn];
                        for (int i = 0; i < nn; i++)
                        {
                            string[] words = lines[i].Split(' ');
                            wxcoor[i].nodeid = words[0].Trim('\t');
                            nlist[i] = wxcoor[i].nodeid;
                            double[] values = new double[3];
                            int count = 0;
                            for (int k = 1; k < words.Length; k++)
                            {
                                if (words[k].Trim().Length == 0) continue;
                                values[count] = Convert.ToDouble(words[k]) * scaleL;
                                count++;
                            }
                            wxcoor[i].x = values[0];
                            wxcoor[i].y = values[1];
                            wxcoor[i].z = values[2];
                        }
                        Debug.WriteLine("Number of Nodes:" + nn.ToString());
                    }
                    else if (node.Name == "section")
                    {
                        if (node.Attributes!.Count > 0) { unitL = node.Attributes[0].Value; }
                        else unitL = "m";
                        scaleL = ToMeters(unitL);
                        string content = node.InnerText.ToString().Trim();
                        string[] lines = content.Split('\n');
                        nsec = lines.Length;
                        wxsection = new Section[nsec];
                        nseclist = new string[nsec];
                        Debug.WriteLine("Number of Sections:" + nsec.ToString());
                        for (int i = 0; i < nsec; i++)
                        {
                            string[] words = lines[i].Split(new char[] { ' ', '=' });
                            string s5 = words[0].Trim('\t');
                            nsection.sectionid = s5;
                            nseclist[i] = s5;
                            for (int k = 0; k < words.Length; k++)
                            {
                                if (words[k] == "Type") nsection.sectiontype = words[k + 1];
                                else if (words[k] == "OD") nsection.od = Convert.ToDouble(words[k + 1]) * scaleL;
                                else if (words[k] == "WTH") nsection.wth = Convert.ToDouble(words[k + 1]) * scaleL;
                                else if (words[k] == "SECID") nsection.sectionid = words[k + 1];
                                if (nsection.sectiontype == "Tube") nsection.SecTube();
                                if (nsection.sectiontype == "AISC")
                                {
                                    scaleL = ToMeters("inch");
                                    scaleF = TokNperM("lbf/ft");
                                    nsection.AISC(scaleL, scaleF);
                                }
                            }
                            wxsection[i] = nsection;
                            Debug.WriteLine("OD" + "\t" + nsection.od.ToString() + "\t" + "WTH" + "\t" + nsection.wth.ToString());
                        }

                    }
                    else if (node.Name == "material")
                    {
                        Material nmaterial = new();

                        if (node.Attributes!.Count > 0) { unitS = node.Attributes[0].Value; unitDen = node.Attributes[1].Value; }
                        else { unitS = "kN/m2"; unitDen = "kN/m3"; }
                        scaleS = TokNperM2(unitS);
                        scaleDen = TokNperM3(unitDen);
                        string content = node.InnerText.ToString().Trim();
                        string[] lines = content.Split('\n');
                        nmat = lines.Length;
                        wxmaterial = new Material[nmat];
                        nmatlist = new string[nmat];
                        for (int i = 0; i < nmat; i++)
                        {
                            string[] words = lines[i].Split(new char[] { ' ', '=' });
                            nmatlist[i] = words[0].Trim('\t');
                            string matname = words[1].Trim('\t');
                            nmaterial.matname = matname;
                            if (matname == "Steel")
                            {
                                nmaterial.matname = "Steel";
                                nmaterial.Ematerial = 200e+06; //kN / m2 or 200GPa;
                                nmaterial.poisson = 0.28;
                                nmaterial.denMaterial = 78.5; //kN / m3;
                                nmaterial.Gmod = 79.3e+06; // kN / m2 or 79.3 GPa;
                            }
                            else if (matname == "General")
                            {
                                for (int k = 0; k < words.Length; k++)
                                {
                                    if (words[k] == "E") nmaterial.Ematerial = Convert.ToDouble(words[k + 1]) * scaleS;
                                    if (words[k] == "Density") nmaterial.denMaterial = Convert.ToDouble(words[k + 1]) * scaleDen;
                                    if (words[k] == "Poisson") Double.TryParse(words[k + 1], out nmaterial.poisson);
                                }
                            }
                            else if (matname == "Titanium")
                            {
                                nmaterial.matname = "Titanium";
                                nmaterial.Ematerial = 113e+06; // !kN / m2 or 120GPa
                                nmaterial.poisson = 0.3;
                                nmaterial.denMaterial = 44.13; // !e + 03 !kN / m3
                                nmaterial.Gmod = 45e+06; // !kN / m2 or 45 GPa

                            }
                            wxmaterial[i] = nmaterial;

                        }
                        Debug.WriteLine("Number of Materials:" + nmat.ToString());
                    }
                    else if (node.Name == "elements")
                    {
                        string content = node.InnerText.ToString().Trim();
                        string[] lines = content.Split('\n');
                        ne = lines.Length;
                        wxelement = new Member[ne];
                        nelemlist = new string[ne];
                        for (int i = 0; i < ne; i++)
                        {
                            string[] words = lines[i].Split(' ');
                            wxelement[i].memberid = words[0].Trim('\t');
                            nelemlist[i] = wxelement[i].memberid;
                            wxelement[i].inc1 = Array.IndexOf(nlist!, (words[1].Trim('\t')));
                            wxelement[i].inc2 = Array.IndexOf(nlist!, (words[2].Trim('\t')));
                            wxelement[i].sec_id = Array.IndexOf(nseclist!, words[3].Trim('\t'));
                            wxelement[i].mat_id = Array.IndexOf(nmatlist!, words[4].Trim(new char[] { '\t', '\r' }));
                            int j1 = wxelement[i].inc1;
                            int j2 = wxelement[i].inc2;
                            wxelement[i].d = wxcoor![j2].distance(wxcoor[j1]);
                            wxelement[i].cx = (wxcoor[j2].x - wxcoor[j1].x) / wxelement[i].d;
                            wxelement[i].cy = (wxcoor[j2].y - wxcoor[j1].y) / wxelement[i].d;
                            wxelement[i].cz = (wxcoor[j2].z - wxcoor[j1].z) / wxelement[i].d;
                            int L = Math.Abs(j2 - j1);
                            if (ms < L) ms = L;

                        }
                        ms = ndf * (ms + 1);
                        tk = new double[n, ms];
                        fem_dload = new double[ne, ndfel];
                        Debug.WriteLine("Number of Elements:" + ne.ToString());
                    }
                    else if (node.Name == "boundary")
                    {
                        string content = node.InnerText.ToString().Trim();
                        string[] lines = content.Split('\n');
                        nbn = lines.Length;
                        wxboundary = new Boundary[nbn];
                        nbnlist = new string[nbn];
                        Debug.WriteLine("Number of Boundary Nodes:" + nbn.ToString());
                        for (int i = 0; i < nbn; i++)
                        {
                            Boundary bnd;// = new();
                            bnd.w = new double[ndf];
                            bnd.ib = new int[ndf + 1];
                            string[] words = lines[i].Split(' ');
                            bnd.nodebid = words[0].Trim('\t');
                            string s1 = words[0].Trim('\t');
                            int n1 = Array.IndexOf(nlist!, s1);
                            bnd.ib[0] = n1;
                            bnd.nbtype = words[1].Trim('\t');
                            bnd.SetBound(structure);
                            wxboundary[i] = bnd;
                            Debug.WriteLine("Node" + "\t" + n1.ToString() + "\t" + bnd.nodebid + "\t" + bnd.nbtype);
                        }
                    }
                    else if (node.Name == "loading")
                    {
                        if (node.Attributes!.Count > 0) { unitL = node.Attributes[0].Value; unitF = node.Attributes[1].Value; }
                        else { unitL = "m"; unitF = "kN"; }
                        scaleL = ToMeters(unitL);
                        scaleF = TokN(unitF);
                        XmlNodeList xcaseList = xDoc.SelectNodes("//case")!;
                        if (xcaseList.Count > 0)
                        {
                            nlc = xcaseList.Count;
                            Debug.WriteLine("Number of Loading Cases:" + xcaseList.Count);
                            for (int klc = 0; klc < xcaseList!.Count; klc++)
                            {
                                if (xcaseList[klc]!.Attributes!.Count > 0)
                                    caseno = xcaseList[klc]!.Attributes![0].Value;
                                Debug.WriteLine("Loading Case:" + caseno);
                                loadtype = xcaseList[klc]!.ChildNodes[0]!.LocalName;
                                string content = (xcaseList[klc]!.InnerText.Trim());
                                string[] lines = content.Split('\n');
                                Debug.WriteLine(loadtype + ": " + lines.Length.ToString());
                                if (loadtype == "loaded-nodes")
                                {
                                    int nldnodes = lines.Length;
                                    for (int i = 0; i < nldnodes; i++)
                                    {
                                        Jload jload = new();
                                        jload.P = new double[ndf];
                                        string[] words = lines[i].Split(new char[] { ' ', '=' });
                                        for (int k = 0; k < words.Length; k++)
                                        {
                                            if (loadtype == "loaded-nodes")
                                            {
                                                if (words[k] == "nodeid")
                                                {
                                                    jload.nodeid = words[k + 1];
                                                    jload.jnodeid = Array.IndexOf(nlist!, jload.nodeid);
                                                    jload.jlc = klc;
                                                    Debug.Write(words[k + 1] + "\t");
                                                }
                                                else if (words[k] == "Px")
                                                {
                                                    jload.P[0] = Convert.ToDouble(words[k + 1]) * scaleF;
                                                    Debug.Write(words[k] + "\t"); Debug.Write(words[k + 1]);
                                                }

                                                else if (words[k] == "Py")
                                                {
                                                    jload.P[1] = Convert.ToDouble(words[k + 1]) * scaleF;
                                                    Debug.Write(words[k] + "\t"); Debug.Write(words[k + 1]);
                                                }

                                                else if (words[k] == "Mz")
                                                {
                                                    jload.P[2] = Convert.ToDouble(words[k + 1]) * scaleF;
                                                    Debug.Write(words[k] + "\t"); Debug.Write(words[k + 1]);
                                                }
                                            }
                                        }
                                        wxjload.Add(jload);
                                    }
                                }
                                if (loadtype == "loaded-members")
                                {
                                    int nldmembers = lines.Length;
                                    for (int i = 0; i < nldmembers; i++)
                                    {
                                        Mload mload = new();
                                        mload.f = new double[ndfel];
                                        mload.vlocal = new double[ndfel];
                                        string[] words = lines[i].Split(new char[] { ' ', '=' });
                                        for (int k = 0; k < words.Length; k++)
                                        {
                                            if (words[k] == "memid")
                                            {
                                                mload.mem_id = words[k + 1];
                                                mload.elemindex = Array.IndexOf(nelemlist!, mload.mem_id);
                                                mload.lc_number = klc;
                                                Debug.WriteLine(words[k + 1]);
                                            }
                                            else if (words[k] == "wa")
                                                mload.wa = Convert.ToDouble(words[k + 1]) * scaleF * scaleL;

                                            else if (words[k] == "wb")
                                                mload.wb = Convert.ToDouble(words[k + 1]) * scaleF * scaleL;

                                            else if (words[k] == "a")
                                                mload.a = Convert.ToDouble(words[k + 1]) * scaleL;
                                            else if (words[k] == "P")
                                                mload.P = Convert.ToDouble(words[k + 1]) * scaleF;
                                            else if (words[k] == "loadtype")
                                                mload.loadtype = words[k + 1];
                                            else if (words[k] == "sysref")
                                            {
                                                if (words[k + 1] == "globx")
                                                { mload!.f[0] = 1; mload!.f[3] = 1; }
                                                else if (words[k + 1] == "globy")
                                                { mload!.f[1] = 1; mload!.f[4] = 1; }
                                                else if (words[k + 1] == "globz")
                                                { mload!.f[2] = 1; mload!.f[5] = 1; }
                                            }
                                        }
                                        wxmload.Add(mload);
                                    }
                                }
                            }
                        }
                    }
                }
                intforc = new double[ne * ndfel * nlc];
                al = new double[n, nlc];
                reac = new double[n, nlc];
                //Matrix al = new(n, nlc);
                //Matrix reac = new(n, nlc);
                //fem_dload = new double[ne, ndfel];
                if (wxjload != null)
                {
                    for (int i = 0; i < wxjload.Count; i++)
                    {
                        Jload jload;
                        jload = (Jload)wxjload![i]!;
                        int nid = jload.jnodeid;
                        int kdsp = ndf * nid;
                        int klc = jload.jlc;
                        al[kdsp, klc] = jload.P![0];
                        al[kdsp + 1, klc] = jload.P![1];
                        al[kdsp + 2, klc] = jload.P![2];
                    }
                    //#if DEBUG
                    //                for (int i = 0; i < n;)
                    //                {
                    //                    for (int j = 0; j < nlc;)
                    //                    {
                    //                        Debug.Write(al[i, j].ToString() + "\t");
                    //                    }
                    //                    Debug.WriteLine("");
                    //                }
                    //#endif
                }
                if (wxmload != null) Mfemgen();
                //treeView1.ExpandAll();
            }
            public double ToMeters(string unitL)
            {
                meters = unitL switch
                {
                    ("inch") => 1.0 / 39.37,
                    ("ft") => 1.0 / 3.28,
                    ("mm") => 0.001,
                    ("cm") => 0.01,
                    ("m") => 1.0,
                    _ => 1,
                };
                return meters;
            }
            public double TokN(string unitF)
            {
                var tokN = unitF switch
                {
                    ("Te") => g,
                    ("kN") => 1.0,
                    ("kips") => 4.4480,
                    ("lbs") => 4448,
                    ("tonf") => 8.8960,
                    ("tonnef") => g,
                    ("N") => 0.001,
                    _ => 1,
                };
                return tokN;
            }
            public double TokNperM(string unitF)
            {
                var tonewtonperm = unitF switch
                {
                    ("kN/m") => 1.0,
                    ("lbf/ft") => 0.014594,
                    ("kip/inch") => 8275.0000,
                    ("kgf/m") => g / 1000,
                    ("tonf/ft") => 29.1900,
                    ("tonnef/m") => g,
                    _ => 1,
                };
                return tonewtonperm;
            }
            public double TokNperM2(string unitF)
            {
                var tonewtonperm2 = unitF switch
                {
                    ("kN/m2") => 1.0,
                    ("lbf/sqf") => 0.04788,
                    ("ksi") => 6894.7600,
                    ("kgf/m2") => g / 1000,
                    ("tonf/sqf") => 95.7600,
                    ("psi") => 6.89476,
                    _ => 1,
                };
                return tonewtonperm2;
            }
            public double TokNperM3(string unitF)
            {
                var tonewtonperm3 = unitF switch
                {
                    ("kN/m3") => 1.0,
                    ("lbf/cuf") => 0.004448,
                    ("kip/inch3") => 271.44714116097,
                    ("kgf/m3") => g / 1000,
                    ("tonf/cuf") => 31.4200,
                    ("tonnef/m3") => g,
                    _ => 1,
                };
                return tonewtonperm3;
            }
            public void Mfemgen()
            {
                double[] vl = new double[ndfel];
                double ra = 0, rb = 0, rma = 0, rmb = 0;
                int nmcases = wxmload.Count;
                for (int i = 0; i < nmcases; i++)
                {
                    Mload mload;
                    mload = (Mload)wxmload[i]!;
                    //mload.f = new double[ndfel];
                    mload.vlocal = new double[ndfel];
                    int mn = mload.elemindex;
                    //Matrix rot = Rotmat(mn);
                    double[,] rot = Rotmat(mn);

                    //int n1 = Array.IndexOf(nlist!, wxelement![mn].inc1);
                    int n1 = wxelement![mn].inc1;
                    //int n2 = Array.IndexOf(nlist!, wxelement[mn].inc2);
                    int n2 = wxelement[mn].inc2;
                    double dl = wxelement[mn].d;

                    //glm::quat pload = { 0, 0,-1, 0 };
                    //glm::quat q = { std::cos(alpha / 2),0,std::sin(alpha / 2) ,0};
                    //glm::quat p1 = q * pload * glm::conjugate(q);

                    int kdsp1 = ndf * n1;
                    int kdsp2 = ndf * n2;
                    int klc = mload.lc_number;
                    double[] f = mload.f!;
                    double a;
                    if (mload.loadtype == "wload")
                    {
                        double wa = mload.wa;
                        double wb = mload.wb;
                        a = mload.a;
                        ra = wa * (Math.Pow((dl - a), 3)) * (dl + a) / (2.0 * dl * dl * dl);
                        rb = ((wa + wb) * (dl - a) / 2.0) - ra;
                        rma = wa * (Math.Pow((dl - a), 3)) * (dl + 3.0 * a) / (12.0 * dl * dl);
                        rmb = (ra * dl) - ((wa * Math.Pow((dl - a), 2)) / 2.0) - ((wb - wa) * (dl - a) * (dl - a) / 6.0) + rma;
                    }
                    else if (mload.loadtype == "pload")
                    {
                        double p = mload.P;
                        a = mload.a;
                        double B = dl * (1 - a);
                        double A = a * dl;
                        ra = p * B * B * (3 * A + B) / (dl * dl * dl);
                        rb = p * A * A * (A + 3 * B) / (dl * dl * dl);
                        rma = p * A * B * B / (dl * dl);
                        rmb = p * A * A * B / (dl * dl);
                    }
                    if (structure == "Frame2D")
                    {
                        vl = MvecMultiply(Transpose(rot), f);
                        //vl = ~rot * f;
                        vl[0] = -ra * vl[0];
                        vl[1] = -ra * vl[1];
                        vl[2] = -rma * vl[2];
                        vl[3] = -rb * vl[3];
                        vl[4] = -rb * vl[4];
                        vl[5] = rmb * vl[5];
                        mload.vlocal = vl;
                    }
                    else if (structure == "Grid")
                    {
                        vl = MvecMultiply(Transpose(rot), f);
                        //vl = ~rot * f;
                        vl[0] = 0.0; vl[2] = -ra * vl[2]; vl[1] = -rma; vl[3] = 0.0; vl[5] = -rb; vl[4] = rmb;
                        mload.vlocal = vl;
                    }
                    else if (structure == "Truss2D")
                    {
                        vl = MvecMultiply(Transpose(rot), f);
                        //vl = ~rot * f;
                        vl[0] = -ra * vl[0]; vl[1] = 0; vl[2] = -rb * vl[2]; vl[3] = 0;
                        mload.vlocal = vl;
                    }
                    else if (structure == "Truss3D")
                    {
                        vl = MvecMultiply(Transpose(rot), f);
                        //vl = ~rot * f;
                        mload.vlocal = vl;

                    }
                    else if (structure == "Frame3D")
                    {
                        vl = MvecMultiply(Transpose(rot), f);
                        //vl = ~rot * f;
                        mload.vlocal = vl;

                    }
                    wxmload[i] = mload;
                    //double[] vglob = ~rot * vl;// MvecMultiply(Transpose(rot), vl);
                    double[] vglob = MvecMultiply(Transpose(rot), vl);
                    for (int j = 0; j < ndf; j++)
                    {
                        al![kdsp1 + j, klc] = al[kdsp1 + j, klc] + vglob[j];
                        al[kdsp2 + j, klc] = al[kdsp2 + j, klc] + vglob[j + ndf];
                    }
                }
            }
            public void forcegen()
            {
                double[] u = new double[ndfel];

                for (int klc = 0; klc < nlc; klc++)
                {
                    for (int nel = 0; nel < ne; nel++)
                    {
                        int n1 = wxelement![nel].inc1;
                        int n2 = wxelement[nel].inc2;
                        //Matrix rot = Rotmat(nel);
                        double[,] rot = Rotmat(nel);
                        int k1 = ndf * n1;
                        int k2 = ndf * n2;
                        int j2;
                        int j1;
                        for (int i = 0; i < ndf; i++)
                        {
                            j1 = k1 + i;
                            j2 = k2 + i;
                            u[i] = al![j1, klc];
                            u[i + ndf] = al[j2, klc];
                        }
                        //double[] ul = rot * u;// MvecMultiply(rot, u);
                        double[] ul = MvecMultiply(rot, u);
                        //Matrix elst = elem_stiff(nel);
                        double[,] elst = elem_stiff(nel);
                        //double[] f = elst * ul;// MvecMultiply(elst, ul);
                        double[] f = MvecMultiply(elst, ul);

                        int i1 = ndfel * nel + ne * ndfel * klc;
                        int i2;
                        if (wxmload != null) // Add contribution from member loads
                        {
                            for (int i = 0; i < wxmload.Count; i++)
                            {
                                Mload mload;
                                mload = (Mload)wxmload[i]!;
                                if (mload.lc_number == klc && mload.elemindex == nel)
                                {
                                    for (int k = 0; k < ndfel; k++)
                                    {
                                        i2 = i1 + k;
                                        f[k] += fem_dload![nel, k] - mload.vlocal![k];
                                        intforc![i2] = f[k];
                                    }
                                }
                                else
                                {
                                    for (int k = 0; k < ndfel; k++)
                                    {
                                        i2 = i1 + k;
                                        f[k] += -fem_dload![nel, k];
                                        intforc![i2] = f[k];
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (int k = 0; k < ndfel; k++)
                            {
                                i2 = i1 + k;
                                f[k] += -fem_dload![nel, k];
                                intforc![i2] = f[k];
                            }

                        }
                        //double[] fg = ~rot * f;// MvecMultiply(Transpose(rot), f);
                        double[] fg = MvecMultiply(Transpose(rot), f);
                        for (int i = 0; i < ndf; i++)
                        {
                            j1 = k1 + i;
                            j2 = k2 + i;
                            reac![j1, klc] = reac[j1, klc] + fg[i];
                            reac[j2, klc] = reac[j2, klc] + fg[i + ndf];
                        }
                    }// for nel
                } // for klc
            }
            double[,] elem_stiff(int nel)
            {
                //Matrix elst = new Matrix(ndfel, ndfel);
                double[,] elst = new double[ndfel, ndfel];
                int imat = wxelement![nel].mat_id;
                int isec = wxelement[nel].sec_id;
                double D = wxelement[nel].d;
                double E = wxmaterial![imat].Ematerial;
                double G = wxmaterial[imat].Gmod;
                double AX = wxsection![isec].A;
                double YZ = wxsection[isec].Ix;
                double YY = wxsection[isec].Iy;
                double YX = wxsection[isec].J;

                if (structure == "Frame2D")
                {
                    elst[0, 0] = E * AX / D;
                    elst[0, 3] = -elst[0, 0];
                    elst[1, 1] = 12.0 * E * YZ / (D * D * D);
                    elst[1, 2] = 6.0 * E * YZ / (D * D);
                    elst[1, 4] = -elst[1, 1];
                    elst[1, 5] = elst[1, 2];
                    elst[2, 1] = elst[1, 2];
                    elst[2, 2] = 4.0 * E * YZ / D;
                    elst[2, 4] = -elst[1, 2];
                    elst[2, 5] = 2.0 * E * YZ / D;
                    elst[3, 0] = elst[0, 3];
                    elst[3, 3] = elst[0, 0];
                    elst[4, 1] = elst[1, 4];
                    elst[4, 2] = elst[2, 4];
                    elst[4, 4] = elst[1, 1];
                    elst[4, 5] = elst[2, 4];
                    elst[5, 1] = elst[1, 5];
                    elst[5, 2] = elst[2, 5];
                    elst[5, 4] = elst[4, 5];
                    elst[5, 5] = elst[2, 2];
                }
                else if (structure == "Frame3D")
                {
                    elst[0, 0] = E * AX / D;
                    elst[6, 0] = -E * AX / D;
                    elst[1, 1] = 12.0 * E * YZ / Math.Pow(D, 3);
                    elst[5, 1] = 6.0 * E * YZ / Math.Pow(D, 2);
                    elst[7, 1] = -12.0 * E * YZ / Math.Pow(D, 3);
                    elst[11, 1] = 6.0 * E * YZ / Math.Pow(D, 2);
                    elst[2, 2] = 12.0 * E * YY / Math.Pow(D, 3);
                    elst[4, 2] = -6.0 * E * YY / Math.Pow(D, 2);
                    elst[8, 2] = -12.0 * E * YY / Math.Pow(D, 3);
                    elst[10, 2] = -6.0 * E * YY / Math.Pow(D, 2);
                    elst[3, 3] = G * YX / D;
                    elst[9, 3] = -G * YX / D;
                    elst[2, 4] = -6.0 * E * YY / Math.Pow(D, 2);
                    elst[4, 4] = 4.0 * E * YY / D;
                    elst[8, 4] = 6.0 * E * YY / Math.Pow(D, 2);
                    elst[10, 4] = 2.0 * E * YY / D;
                    elst[1, 5] = 6.0 * E * YZ / Math.Pow(D, 2);
                    elst[5, 5] = 4.0 * E * YZ / D;
                    elst[7, 5] = -6.0 * E * YZ / Math.Pow(D, 2);
                    elst[11, 5] = 2.0 * E * YZ / D;
                    elst[0, 6] = -E * AX / D;
                    elst[6, 6] = E * AX / D;
                    elst[1, 7] = -12.0 * E * YZ / Math.Pow(D, 3);
                    elst[5, 7] = -6.0 * E * YZ / Math.Pow(D, 2);
                    elst[7, 7] = 12.0 * E * YZ / Math.Pow(D, 3);
                    elst[11, 7] = -6.0 * E * YZ / Math.Pow(D, 2);
                    elst[2, 8] = -12.0 * E * YY / Math.Pow(D, 3);
                    elst[4, 8] = 6.0 * E * YY / Math.Pow(D, 2);
                    elst[8, 8] = 12.0 * E * YY / Math.Pow(D, 3);
                    elst[10, 8] = 6.0 * E * YY / Math.Pow(D, 2);
                    elst[3, 9] = -G * YX / D;
                    elst[9, 9] = G * YX / D;
                    elst[2, 10] = -6.0 * E * YY / Math.Pow(D, 2);
                    elst[4, 10] = 2.0 * E * YY / D;
                    elst[8, 10] = 6.0 * E * YY / Math.Pow(D, 2);
                    elst[10, 10] = 4.0 * E * YY / D;
                    elst[1, 11] = 6.0 * E * YZ / Math.Pow(D, 2);
                    elst[5, 11] = 2.0 * E * YZ / D;
                    elst[7, 11] = -6.0 * E * YZ / Math.Pow(D, 2);
                    elst[11, 11] = 4.0 * E * YZ / D;
                }
                else if (structure == "Truss2D")
                {
                    elst[0, 0] = E * AX / D;
                    elst[1, 0] = 0;
                    elst[0, 1] = -E * AX / D;
                    elst[2, 2] = 0;
                }
                else if (structure == "Truss3D")
                {
                    elst[0, 0] = E * AX / D;
                    elst[3, 0] = -E * AX / D;
                    elst[0, 3] = -E * AX / D;
                    elst[3, 3] = E * AX / D;
                }
                else if (structure == "Grid")
                {
                    elst[0, 0] = G * YX / D;
                    elst[0, 3] = -G * YX / D;
                    elst[1, 1] = 4.0 * E * YZ / D;
                    elst[1, 2] = -6.0 * E * YZ / (D * D);
                    elst[1, 4] = 2.0 * E * YZ / D;
                    elst[1, 5] = 6.0 * E * YZ / (D * D);
                    elst[2, 1] = elst[1, 2];
                    elst[2, 2] = 12.0 * E * YZ / (D * D * D);
                    elst[2, 4] = elst[1, 2];
                    elst[2, 5] = -12.0 * E * YZ / (D * D * D);
                    elst[3, 0] = -elst[0, 0];
                    elst[3, 3] = elst[0, 0];
                    elst[4, 1] = 2.0 * E * YZ / D;
                    elst[4, 2] = -6.0 * E * YZ / (D * D);
                    elst[4, 4] = 4.0 * E * YZ / D;
                    elst[4, 5] = 6.0 * E * YZ / (D * D);
                    elst[5, 1] = 6.0 * E * YZ / (D * D);
                    elst[5, 2] = -12.0 * E * YZ / (D * D * D);
                    elst[5, 4] = 6.0 * E * YZ / (D * D);
                    elst[5, 5] = 12.0 * E * YZ / (D * D * D);
                }
                return elst;
            }
            public void elassgen(int nel, double[,] elst)
            {
                int N1 = wxelement![nel].inc1;// con[0];// inc1;// ele[NEL - 1].n1;
                int N2 = wxelement![nel].inc2;// con[1];// inc2;// ele[NEL - 1].n2;
                int J1 = ndf * N1;
                int J2 = ndf * N2;
                int KR;
                int KC;
                /*------------------------------------
                * NEL	=number of the current node
                * N1	=number of the start node
                * N2	=number of the END node
                ------------------------------------*/
                try
                {
                    for (int I = 0; I < ndf; I++)//'N1'
                    {
                        KR = J1 + I;
                        for (int J = I; J < ndf; J++)
                        {
                            KC = J - I;
                            tk![KR, KC] = tk[KR, KC] + elst[I, J];
                        }
                    }

                    for (int I = 0; I < ndf; I++)//'N2'
                    {
                        KR = J2 + I;
                        for (int J = I; J < ndf; J++)
                        {
                            KC = J - I;
                            tk![KR, KC] = tk[KR, KC] + elst[I + ndf, J + ndf];
                        }
                    }

                    int IC;
                    if (N1 < N2)
                    {
                        for (int I = 0; I < ndf; I++)//'N1<N2'
                        {
                            KR = J1 + I;
                            IC = J2 - KR;
                            for (int J = 0; J < ndf; J++)
                            {
                                KC = IC + J;
                                tk![KR, KC] = tk[KR, KC] + elst[I, J + ndf];
                            }
                        }
                    }
                    if (N1 > N2)
                    {
                        for (int I = 0; I < ndf; I++)
                        {
                            KR = J2 + I;
                            IC = J1 - KR;
                            for (int J = 0; J < ndf; J++)
                            {
                                KC = IC + J;
                                tk![KR, KC] = tk[KR, KC] + elst[I + ndf, J];
                            }
                        }
                    }
                }
                catch (Exception ex)
                {
                    Debug.WriteLine(ex.Message.ToString());
                    Debug.WriteLine("Program Exception " + N1.ToString() + "-" + N2.ToString());
                    Debug.WriteLine("Element " + nel.ToString());
                    //Debug.WriteLine("KR=" + KR.ToString());
                    //Debug.WriteLine("KC=" + KC.ToString());
                    Debug.WriteLine(ex.ToString());
                }



            }
            public void bound()
            {
                // introduction of the boundary conditions
                for (int l = 0; l < nbn; l++) //do 100 l=1,int
                {
                    int no = wxboundary![l].ib![0];
                    int k1 = ndf * no;
                    for (int i = 0; i < ndf; i++)
                    {
                        if (wxboundary[l].ib![i] == 0)// set diagonal coefficient of tk equal to 1  and place prescribed value in al
                        {
                            int kr = k1 + i;
                            tk![kr, 0] = 1.0;
                            for (int icol = 0; icol < nlc; icol++)
                                al![kr, icol] = reac![kr, icol];
                            for (int j = 1; j < ms; j++) //modify values in columns
                            {
                                int kv = kr + j - 1;
                                if (kv < n - 1)  // modify row of tk and corresponding elements in al
                                {
                                    for (int icol = 0; icol < nlc; icol++)
                                    {
                                        al![kv, icol] = al[kv, icol] - tk[kr, j] * reac![kr, icol];
                                    }
                                    tk[kr, j] = 0.0;
                                }
                                if (kr > j)   // modify column in tk and corresponding element in al
                                {
                                    kv = kr - j + 1;
                                    for (int icol = 0; icol < nlc; icol++)
                                        al![kv, icol] = al[kv, icol] - tk[kv, j] * reac![kr, icol];
                                    tk[kv, j] = 0.0;
                                }
                            } //for j
                        }// if ib
                    } // for i
                }// for l
#if DEBUG
                using StreamWriter Debugfile = new(filelog, append: true);
                Debugfile.WriteLine("Stiffness Matrix after Bound");
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < ms; j++)
                    {
                        Debugfile.Write(tk![i, j].ToString() + "\t");
                    }
                    Debugfile.WriteLine("");
                }
#endif
            }
            public void bound2()
            {
                ///introduction of the boundary conditions
                for (int l = 1; l <= nbn; l++) //do 100 l=1,nbn
                {
                    ///no=number of the current boundary

                    int no = wxboundary![l - 1].ib![0] + 1;
                    int k1 = ndf * (no - 1);
                    int i;
                    for (i = 1; i <= ndf; i++) //do 100 i=1,ndf
                    {
                        //IF(IB(L2)) 100,10,100
                        if (wxboundary[l - 1].ib![i] is < 0 or > 0) goto L100;
                        ///set diagonal coefficient of tk equal to 1 and place prescribed value in al
                        int kr = k1 + i;
                        int icol;
                        for (int j = 2; j <= ms; j++) //do 50 j=2,ms
                        {
                            int kv = kr + j - 1;
                            //if(n-kv)30,20,20
                            if (n < kv) goto l30;
                            ///modify row of tk and corresponding elements in al
                            //20; //DO 60 ICOL=1,NLC
                            for (icol = 1; icol <= nlc; icol++)
                                al![kv - 1, icol - 1] = al[kv - 1, icol - 1] - tk![kr - 1, j - 1] * reac![kr - 1, icol - 1];
                            tk![kr - 1, j - 1] = 0.0;
                        l30: kv = kr - j + 1;
                            //IF(KV)50,50,40
                            if (kv <= 0) goto l50;
                            ///modify column in tk and corresponding element in al
                            //40; //DO 70 ICOL=1,NLC
                            for (icol = 1; icol <= nlc; icol++)
                                al![kv - 1, icol - 1] = al[kv - 1, icol - 1] - tk![kv - 1, j - 1] * reac![kr - 1, icol - 1];
                            tk![kv - 1, j - 1] = 0.0;
                        //l50: //continue do
                        l50: tk![kr - 1, 0] = 1.0;
                        }
                        for (icol = 1; icol <= nlc; icol++) //do 200 icol=1,nlc
                            al![kr - 1, icol - 1] = reac![kr - 1, icol - 1];
                        L100:; //Continue Do
                    }
                }
                return;
            }
            public void OUTPT2(int kiter, string title, string ptitle)
            {
                /// <summary>
                ///PROGRAM TO OUTPUT JOINT DISPLACEMENT, NODAL REACTION AND MEMBER ForceS
                ///</summary>
                string LINE;
                double[] Ulen = new double[ndf + 1], Ufor = new double[ndf + 1], rmom = new double[ndfel + 1];
                double sLen, sFor;

                int I, J, K, K1, K2, KLC, L1, NO, NEL, N1, count;

                string Tab;
                string dispheader = string.Format("\t\t" + "Node\t" + "dx\t" + "dy\t" + "rotz");
                string dispheader2 = string.Format("\t\t" + "(" + unitL + ")" + "(" + unitL + ")" + "(rad)");
                string forheader = string.Format("\t\t" + "Node\t" + "Px\t" + "Py\t" + "Mz");
                string forheader2 = string.Format("{0,10} {1,10} {2,10} {3,10}", "    ", "<" + unitF + ">", "<" + unitF + ">", "<" + unitF + unitL + ">");
                string memheader = string.Format("{0,10} {1,10} {2,10} {3,10} {4,10}", "Element", "Node", "Px", "Py", "Mz");
                string memheader2 = string.Format("{0,10} {1,10} {2,10} {3,10} {4,10}", "  ", "    ", "<" + unitF + ">", "<" + unitF + ">", "<" + unitF + unitL + ">");
                Tab = "       ";
                // --------------------------------------------------------------------
                // WRITE NODAL DISPLACEMENTS
                // ---------------------------------------------------------------------
                using StreamWriter file = new(fileout);

                //KITER = 1;
                file.WriteLine("\t\t" + ptitle);
                LINE = title!;
                file.WriteLine(title);

                DateTime now = DateTime.Now;
                file.WriteLine(now);
                file.WriteLine("-----------------------------------------------------");
                if ((kiter > 0))
                    file.WriteLine("Number of Iterations " + kiter.ToString());
                sLen = 1.0;// ToMeters(unitL!);
                sFor = 1;// TokN(unitF!);
                         //fyield = 36.0 * TokNperM2("ksi");
                Ulen[0] = 1.0 / sLen;
                Ulen[1] = 1.0 / sLen;
                Ulen[2] = 1.0;
                Ufor[0] = 1.0 / sFor;
                Ufor[1] = 1.0 / sFor;
                Ufor[2] = 1.0 / (sFor * sLen);
                for (KLC = 1; KLC <= nlc; KLC++)
                {
                    file.WriteLine("Nodal Displacements for Loading " + KLC.ToString());
                    file.WriteLine("Active Units: " + unitL);
                    if (structure == "3DFrame")
                        file.WriteLine("NODE" + Tab + "DX" + Tab + "DY" + Tab + "DZ" + "ROTX" + "ROTY" + "ROTZ");
                    else if (structure == "2DFrame")
                    {
                        file.WriteLine(dispheader);
                        file.WriteLine(dispheader2);
                    }
                    for (I = 1; I <= nn; I++)
                    {
                        K1 = ndf * (I - 1) + 1;
                        K2 = K1 + ndf - 1;
                        file.Write("{0,12}\t", nlist![I - 1]);
                        count = 0;
                        for (J = K1; J <= K2; J++)
                        {
                            file.Write("{0,12}\t", String.Format("{0:0.0000}", al![J - 1, KLC - 1] * Ulen[count]));
                            count += 1;
                        }
                        file.WriteLine("");
                    }
                }
                // ------------------------------------------------------
                // WRITE NODAL REACTIONS
                //----------------------------------------------------
                file.WriteLine("-----------------------------------------------------");
                for (KLC = 1; KLC <= nlc; KLC++)
                {
                    file.WriteLine("Nodal Reactions for Loading " + KLC.ToString());
                    file.WriteLine("Active Units: " + unitF + " " + unitL);
                    if (structure == "3DFrame")
                        file.WriteLine("NODE" + Tab + "PX" + Tab + "PY" + Tab + "PZ" + "MX" + Tab + "MY" + Tab + "MZ");
                    else if (structure == "2DFrame")
                    {
                        file.WriteLine(forheader);
                        file.WriteLine(forheader2);
                    }
                    for (I = 1; I <= nbn; I++)
                    {
                        L1 = (ndf + 1) * (I - 1) + 1;
                        NO = wxboundary![I - 1].ib![0];
                        //K1 = Ndf * (NO - 1) + 1;
                        K1 = ndf * NO;// (NO - 1) + 1;
                        K2 = K1 + ndf - 1;
                        file.Write("{0,10}", nlist![NO]);// - 1]);
                        count = 0;
                        for (J = K1; J <= K2; J++)
                        {
                            file.Write("{0,12}", String.Format("{0:#.00}", reac![J, KLC - 1] * Ufor[count]));
                            count += 1;
                        }
                        file.WriteLine("");
                    }
                }
                file.WriteLine("-----------------------------------------------------");
                //---------------------------------------------------
                // OUTPUT MEMBER END ForceS
                // ---------------------------------------------------
                //103:    'IF(KGO.LT.4) RETURN
                for (KLC = 1; KLC <= nlc; KLC++)
                {
                    file.WriteLine("Member End Forces for Loading " + KLC.ToString());
                    file.WriteLine("Active Units: " + unitL);
                    if (structure == "3DFrame")
                        file.WriteLine("ELEMENT" + Tab + "NODE     " + "FX" + Tab + "FY" + Tab + "FZ" + Tab + "MX" + Tab + "MY" + Tab + "MZ");
                    else if (structure == "2DFrame")
                    {
                        file.WriteLine(memheader);
                        file.WriteLine(memheader2);
                    }
                    for (NEL = 1; NEL <= ne; NEL++)
                    {
                        K1 = ndfel * (NEL - 1) + ne * ndfel * (KLC - 1);
                        K2 = K1 + 2;
                        N1 = nne * (NEL - 1);
                        file.Write("{0,10}", String.Format(wxelement![NEL - 1].memberid!));
                        file.Write("{0,10}", String.Format(wxelement[NEL - 1].inc1.ToString()));//[CON[N1]]));
                        count = 0;
                        for (K = K1; K <= K2; K++)
                        {
                            file.Write("{0,12}", String.Format("{0:#.00}", intforc![K] * Ufor[count]));
                            rmom[count] = intforc[K];
                            count += 1;
                        }
                        // rmom(0) = FORC(K1)
                        // rmom(1) = FORC(K2)
                        file.WriteLine("");
                        LINE = NEL.ToString();
                        K1 = K2 + 1;
                        K2 = K1 + 2;
                        file.Write("{0,10}", String.Format(wxelement[NEL - 1].memberid));
                        file.Write("{0,10}", String.Format(wxelement[NEL - 1].inc2.ToString()));//[CON[N1 + 1]]));
                        count = 0;
                        for (K = K1; K <= K2; K++)
                        {
                            file.Write("{0,12}", String.Format("{0:#.00}", intforc![K] * Ufor[count]));
                            rmom[count + ndf] = intforc[K];
                            count += 1;
                        }
                        // rmom(3) = FORC(K1)
                        // rmom(4) = FORC(K2)
                        file.WriteLine("");
                        //if (code == "AISC")
                        //{
                        //    file.WriteLine("");
                        //    unityratio = 1.0;// AISC(NEL, PROP, rmom, fyield, ref eqt);
                        //    file.WriteLine("{0,10} {1,10} {2,10} {3,12} {4,10}", wxelement![NEL - 1], code, eqt, "URatio", String.Format("{0:0.000}", unityratio));
                        //    file.WriteLine("");
                        //}
                        //else if (code == "NDS")
                        //{
                        //    file.WriteLine("");
                        //    unityratio = 1.0;// NDS(NEL, PROP, rmom, ref eqt);
                        //    file.WriteLine("{0,10} {1,10} {2,10} {3,12} {4,10}", wxelement![NEL - 1], code, eqt, "URatio", String.Format("{0:0.000}", unityratio));
                        //    file.WriteLine("");
                        //}

                    }
                    file.WriteLine("-----------------------------------------------------");
                }
                DateTime TimeStr = DateTime.Now;
                file.WriteLine("{0:D}", TimeStr);
                file.WriteLine();
                //file.Close();
            }
            public static double[,] BTAB3(double[,] A, double[,] B)
            {
                int N, I, J, K;
                N = A.GetUpperBound(0) + 1;
                double[] V = new double[N];
                /// C   THIS PROGRAM COMPUTES THE MATRIX OPERATION
                /// A = TRANSPOSE(B) * A * B and store in A
                /// N = ACTUAL ORDER Of A And B
                /// V = AUXILIARY VECTOR

                for (I = 1; I <= N; I++) // DO 10 I=1,N
                {
                    for (J = 1; J <= N; J++) // DO 5 J=1,N
                    {
                        V[J - 1] = 0.0;
                        for (K = 1; K <= N; K++) // DO 5 K=1,N
                            V[J - 1] = V[J - 1] + A[I - 1, K - 1] * B[K - 1, J - 1];
                    }
                    for (J = 1; J <= N; J++) // DO 10 J=1,N
                        A[I - 1, J - 1] = V[J - 1];
                }
                /// COMPUTE TRANSPOSE(B) * A And STORE IN A
                for (J = 1; J <= N; J++) // DO 20 J=1,N
                {
                    for (I = 1; I <= N; I++) // DO 15 I=1,N
                    {
                        V[I - 1] = 0.0;
                        for (K = 1; K <= N; K++) // DO 15 K=1,N
                            V[I - 1] = V[I - 1] + B[K - 1, I - 1] * A[K - 1, J - 1];
                    }
                    for (I = 1; I <= N; I++) // DO 20 I=1,N
                        A[I - 1, J - 1] = V[I - 1];
                }
                return A;
            }
            public static double[] MvecMultiply(double[,] mat, double[] vec)
            {
                int n = mat.GetUpperBound(0) + 1;
                int m = mat.GetUpperBound(1) + 1;
                double[] v = new double[n];
                //Parallel.For(0, n, i =>
                for (int i = 0; i < n; i++)
                {
                    v[i] = 0.0;
                    for (int j = 0; j < m; j++)
                    {
                        v[i] += mat[i, j] * vec[j];
                    }
                }//);
                return v;
            }
            public static double[,] Transpose(double[,] mat)
            {
                int n = mat.GetUpperBound(0) + 1;
                int m = mat.GetUpperBound(1) + 1;
                double[,] v = new double[n, m];
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < m; j++)
                    {
                        v[i, j] = mat[j, i];
                    }
                }
                return v;
            }
            public void dloadgen()
            {
                int isec, imat, inc1, inc2, kdsp1, kdsp2;
                double[] f = new double[ndfel]; double[] vlocal = new double[ndfel];
                double wload, ra, rma, dl, totalwht;
                totalwht = 0.0;
                for (int nel = 1; nel <= ne; nel++)
                {
                    double[,] rot = new double[ndfel, ndfel];
                    isec = wxelement![nel].sec_id;
                    imat = wxelement[nel].mat_id;
                    inc1 = wxelement[nel].inc1;
                    inc2 = wxelement[nel].inc2;
                    dl = wxelement[nel].d;
                    wload = wxsection![isec].A * wxmaterial![imat].denMaterial;
                    kdsp1 = ndf * inc1;
                    kdsp2 = ndf * inc2;
                    ra = wload * dl / 2.0;
                    rma = wload * dl * dl / 12.0;
                    totalwht += wload * dl;
                    if (structure == "Frame3D")
                    {
                        f[1] = 0.0; f[2] = -1.0; f[3] = 0.0; f[4] = 0.0; f[5] = 0.0; f[6] = 0.0;
                        f[7] = 0.0; f[8] = -1.0; f[9] = 0.0; f[10] = 0.0; f[11] = 0.0; f[12] = 0.0;
                        rot = Rotmat(nel);
                        vlocal = MvecMultiply(rot, f);
                        vlocal[1] = ra * vlocal[1];
                        vlocal[2] = ra * vlocal[2];
                        vlocal[3] = ra * vlocal[3];
                        vlocal[4] = -rma * vlocal[4];
                        vlocal[5] = -rma * vlocal[5];
                        vlocal[6] = -rma * vlocal[6];
                        vlocal[7] = ra * vlocal[7];
                        vlocal[8] = ra * vlocal[8];
                        vlocal[9] = ra * vlocal[9];
                        vlocal[10] = rma * vlocal[10];
                        vlocal[11] = rma * vlocal[11];
                        vlocal[12] = rma * vlocal[12];
                    }
                    else if (structure == "Frame2D")
                    {
                        f[1] = 0.0; f[2] = -1.0; f[3] = 0.0; f[4] = 0.0; f[5] = -1.00; f[6] = 0.0;// !dload unit Y vector
                        rot = Rotmat(nel);
                        vlocal = MvecMultiply(rot, f);
                        vlocal[1] = ra * vlocal[1];
                        vlocal[2] = ra * vlocal[2];
                        vlocal[3] = -rma * vlocal[3];
                        vlocal[4] = ra * vlocal[4];
                        vlocal[5] = ra * vlocal[5];
                        vlocal[6] = rma * vlocal[6];
                    }
                    else if (structure == "Truss2D")
                    {
                        f[1] = 0.0; f[2] = -1.0; f[3] = 0.0; f[4] = -1.0;
                        rot = Rotmat(nel);
                        vlocal = MvecMultiply(rot, f);
                        vlocal[1] = ra * vlocal[1];
                        vlocal[2] = ra * vlocal[2];
                        vlocal[3] = ra * vlocal[3];
                        vlocal[4] = ra * vlocal[4];
                    }
                    else if (structure == "Truss3D")
                    {
                        f[1] = 0.0; f[2] = -1.0; f[3] = 0.0; f[4] = 0.0; f[5] = -1.0; f[6] = 0.0;
                        rot = Rotmat(nel);
                        vlocal = MvecMultiply(rot, f);
                        vlocal[1] = ra * vlocal[1];
                        vlocal[2] = ra * vlocal[2];
                        vlocal[3] = ra * vlocal[3];
                        vlocal[4] = ra * vlocal[4];
                        vlocal[5] = ra * vlocal[5];
                        vlocal[6] = ra * vlocal[6];
                    }
                    double[] vglob = MvecMultiply(Transpose(rot), vlocal);
                    for (int j = 0; j < ndfel; j++)
                        fem_dload![nel, j] = fem_dload[nel, j] - vlocal[j];

                    for (int klc = 0; klc < nlc; klc++)
                    {
                        for (int j = 0; j < ndf; j++)
                        {
                            al![kdsp1 + j, klc] = al[kdsp1 + j, klc] + vglob[j];
                            al[kdsp2 + j, klc] = al[kdsp2 + j, klc] + vglob[j + ndf];
                        }
                    }
                }
            }
            public void AISC_360_16_ASD()
            {
                double r = 0, area, rmomi, rmomj, rmom = 0.0, sx, sy, rx, ry, Pr = 0.0, slendy = 0.0, fb, cc,
                    Fe, uratio, blngth, E, cc1, Fcr, Pn, alpha = 1.6, zx, zy, Mnx = 0.0, Mrx = 0.0;
                int klc, n1, n2, k1, k2;
                using StreamWriter outfile = new(fileout);

                outfile.WriteLine("ANSI/AISC Code 360-16 (ASD/LRFD) - July 7, 2016 - Revised June 2019 ");

                for (klc = 1; klc <= nlc; klc++)
                {
                    outfile.WriteLine("Load Case: " + klc.ToString());
                    outfile.WriteLine("Element\t" + "Unity Ratio\t" + "Equation");
                    for (int nel = 1; nel <= ne; nel++)
                    {
                        //strcpy(label,elemid[NEL - 1].c_str());
                        n1 = wxelement![nel - 1].inc1;
                        n2 = wxelement[nel - 1].inc2;
                        k1 = ndfel * (nel - 1) + ne * ndfel * (klc - 1);
                        k2 = k1 + 2;
                        blngth = wxelement[nel - 1].d;
                        int imat = wxelement[nel - 1].mat_id;
                        E = wxmaterial![imat].Ematerial;
                        int isec = wxelement[nel - 1].sec_id;
                        area = wxsection![isec].A;
                        sx = wxsection[isec].Sx;
                        sy = wxsection[isec].Sy;
                        rx = wxsection[isec].rx;
                        ry = wxsection[isec].ry;
                        zx = wxsection[isec].Zx; zy = wxsection[isec].Zy;
                        double fx;
                        double fy;
                        double faw;
                        if (structure == "Frame2D")
                        {
                            Pr = intforc![k1];
                            rmomi = intforc[k2];
                            k1 = k2 + 1;
                            k2 = k1 + 2;
                            rmomj = intforc[k2];
                            r = Math.Min(rx, ry);
                            rmom = Math.Max(Math.Abs(rmomi), Math.Abs(rmomj));
                            Mnx = fyield * zx / alpha;
                            Mrx = rmom;
                            fx = rmom / sx;             // Compute working stresses
                            faw = Pr / area;
                            slendy = blngth / r;
                            fy = 0.0;
                        }
                        if (structure == "Frame3D")
                        {
                            //Compute working stresses
                            //-----
                            double Mny = fyield * zy / alpha;
                            fx = rmom / sx;
                            fy = rmom / sy;
                            faw = Pr / area;
                            slendy = blngth / r;
                        }
                        //---------------------------
                        // Compute allowable stresses
                        //*----------------------------/
                        fb = 0.66 * fyield;
                        cc = 4.71 * Math.Sqrt(E / fyield);
                        Fe = 5.149359 * E / (Math.Pow(slendy, 2));
                        cc1 = fyield / Fe;
                        if ((slendy <= cc || cc1 <= 2.25)) Fcr = fyield * Math.Pow(0.658, cc1);
                        else Fcr = 0.877 * Fe;
                        Pn = Fcr * area / alpha;
                        if ((Pr < 0.0))
                        {
                            if ((Math.Abs(alpha * Pr / Pn) < 0.2))
                            {
                                uratio = Math.Abs(Pr / (2 * Pn)) + Math.Abs(Mrx / Mnx);// Compute Unity Ratio
                            }
                            else
                            {
                                uratio = Math.Abs(Pr / Pn) + Math.Abs(Mrx / Mnx) * 8 / 9;// Compute Unity Ratio
                            }
                        }
                        else
                        {
                            uratio = Math.Abs(Pr / Pn) + Math.Abs(Mrx / Mnx);// Compute Unity Ratio
                        }
                        outfile.WriteLine(wxelement[nel - 1].memberid.ToString() + "\t\t" + "uratio" + "\t\t" + "eqtn");
                    }
                    outfile.WriteLine("");
                }
            }// AISC
        }
    }
}
