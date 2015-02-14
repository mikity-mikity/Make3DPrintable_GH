using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using mikity.GeometryProcessing;
using ShoNS.Array;
using Rhino.Geometry;
namespace mikity.ghComponents
{
    public partial class MeshTest : Grasshopper.Kernel.GH_Component
    {
        public Mesh inMesh;
        public Mesh outMesh;
        public int select = 0;
        public List<List<Line>> listNormal;
        private DoubleArray computeGradient(int numvar, MeshStructure MS, List<Point3d> nodes, DoubleArray x,Vector3d[] normalPerFaces)
        {
            DoubleArray grad = DoubleArray.Zeros(1, numvar);
            foreach (var v in MS.vertices)
            {
                foreach (var e in v.onering)
                {
                    var N = normalPerFaces[e.owner.N];
                    int index=MS.vertices.IndexOf(v);
                    var cx=x[index*3+0,0];
                    var cy=x[index*3+1,0];
                    var cz=x[index*3+2,0];
                    grad[0, index * 3 + 0] += 2 * N.X * N.X * cx + 2 * N.X * N.Y * cy + 2 * N.X * N.Z * cz - 2 * N.X;
                    grad[0, index * 3 + 1] += 2 * N.Y * N.Y * cy + 2 * N.Y * N.Z * cz + 2 * N.Y * N.X * cx - 2 * N.Y;
                    grad[0, index * 3 + 2] += 2 * N.Z * N.Z * cz + 2 * N.Z * N.X * cx + 2 * N.Z * N.Y * cy - 2 * N.Z;
                }
            }
            return grad;
        }
        private DoubleArray computeResidual(int numvar, int numcon, MeshStructure MS, List<Point3d> nodes, DoubleArray x)
        {
            DoubleArray res = DoubleArray.Zeros(numcon, 1);
            int conoffset = 0;
            foreach (var v in MS.vertices)
            {
                int index = MS.vertices.IndexOf(v);
                //det=0
                foreach (var e in v.star)
                {
                    if (e.P.N >= e.next.P.N) continue;
                    var Q = nodes[e.next.P.N];
                    var P = nodes[e.P.N];
                    var ax = (P - Q).X;
                    var ay = (P - Q).Y;
                    var az = (P - Q).Z;
                    int indexB = MS.vertices.IndexOf(e.next.P);
                    int indexC = MS.vertices.IndexOf(v);
                    var bx = x[indexB * 3 + 0];
                    var by = x[indexB * 3 + 1];
                    var bz = x[indexB * 3 + 2];
                    var cx = x[indexC * 3 + 0];
                    var cy = x[indexC * 3 + 1];
                    var cz = x[indexC * 3 + 2];
                    res[conoffset, 0] = bx * (ay * cz - az * cy) + by * (az * cx - ax * cz) + bz * (ax * cy - ay * cx);
                    conoffset++;
                }
            }
            return res;
        }
        private DoubleArray computeJacob(int numvar, int numcon, MeshStructure MS, List<Point3d> nodes, DoubleArray x)
        {
            DoubleArray jacob = DoubleArray.Zeros(numcon, numvar);
            int conoffset = 0;
            foreach (var v in MS.vertices)
            {
                int index = MS.vertices.IndexOf(v);
                //det=0
                foreach (var e in v.star)
                {

                    if (e.P.N >= e.next.P.N) continue;
                    var Q = nodes[e.next.P.N];
                    var P = nodes[e.P.N];
                    var ax = (P - Q).X;
                    var ay = (P - Q).Y;
                    var az = (P - Q).Z;
                    int indexB = MS.vertices.IndexOf(e.next.P);
                    int indexC = MS.vertices.IndexOf(v);
                    var bx = x[indexB * 3 + 0];
                    var by = x[indexB * 3 + 1];
                    var bz = x[indexB * 3 + 2];
                    var cx = x[indexC * 3 + 0];
                    var cy = x[indexC * 3 + 1];
                    var cz = x[indexC * 3 + 2];
                    jacob[conoffset, indexB * 3 + 0] = (ay * cz - az * cy);
                    jacob[conoffset, indexB * 3 + 1] = (az * cx - ax * cz);
                    jacob[conoffset, indexB * 3 + 2] = (ax * cy - ay * cx);
                    jacob[conoffset, indexC * 3 + 0] = by * az - bz * ay;
                    jacob[conoffset, indexC * 3 + 1] = bz * ax - bx * az;
                    jacob[conoffset, indexC * 3 + 2] = bx * ay - by * ax;
                    conoffset++;
                }
            }
            return jacob;
        }
        public MeshTest()
            : base("MeshTest", "MeshTest", "MeshTest", "Kapybara3D", "Computation")
        {
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("e41842d3-325e-4ed4-9f42-9d125f41e555"); }
        }
        protected override void RegisterInputParams(Grasshopper.Kernel.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("mesh", "M", "mesh", Grasshopper.Kernel.GH_ParamAccess.item);
            pManager.AddIntegerParameter("N", "N", "N", Grasshopper.Kernel.GH_ParamAccess.item);
        }
        protected override void RegisterOutputParams(Grasshopper.Kernel.GH_Component.GH_OutputParamManager pManager)
        {
        }
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            Mesh myMesh=new Mesh();
            if (!DA.GetData(0, ref myMesh)) { return; }
            if (!DA.GetData(1, ref select)) { return; }
            mikity.GeometryProcessing.MeshStructure MS = MeshStructure.CreateFrom(myMesh);
            int numvar = MS.vertices.Count * 3;
            int numcon = 0;
            foreach (var v in MS.vertices)
            {
                foreach (var e in v.star)
                {
                    if (e.P.N < e.next.P.N)
                        numcon++;
                }
            }
            DoubleArray x = DoubleArray.Zeros(numvar, 1);
            List<Point3d> nodes = new List<Point3d>();
            Vector3d[] normalPerFaces = new Vector3d[MS.nFaces];
            nodes.AddRange(myMesh.Vertices.ToPoint3dArray());
            listNormal = new List<List<Line>>();

            //initial guess
            var lN=new List<Line>();
            listNormal.Add(lN);
            foreach (var v in MS.vertices)
            {
                int index = MS.vertices.IndexOf(v);
                Vector3d N = new Vector3d();
                foreach (var e in v.onering)
                {
                    //compute normal
                    var P = nodes[e.P.N];
                    var Q = nodes[e.next.P.N];
                    var R = nodes[e.next.next.P.N];
                    var b = P - Q;
                    var c = R - Q;
                    var n = Vector3d.CrossProduct(b, c);
                    n.Unitize();
                    normalPerFaces[e.owner.N] = n;
                    N += n;
                }
                x[index * 3 + 0, 0] = N[0];
                x[index * 3 + 1, 0] = N[1];
                x[index * 3 + 2, 0] = N[2];
                lN.Add(new Line(nodes[v.N], nodes[v.N] + N));
            }
            for (int i = 0; i < 10; i++)
            {
                var grad = computeGradient(numvar, MS, nodes, x, normalPerFaces);
                var jacob = computeJacob(numvar, numcon, MS, nodes, x);
                var S = jacob.Multiply(jacob.T);
                var solve = new ShoNS.Array.SVD(S.T);
                //if (solve.Rank() < numcon) break;
                
                var lambda = solve.Solve(jacob.Multiply(-grad.T)).T;
                var projGrad = grad + lambda.Multiply(jacob);
                foreach (var v in MS.vertices)
                {
                    int index = MS.vertices.IndexOf(v);
                    x[index * 3 + 0, 0] -= projGrad[0, index * 3 + 0]*0.05;
                    x[index * 3 + 1, 0] -= projGrad[0, index * 3 + 1]*0.05;
                    x[index * 3 + 2, 0] -= projGrad[0, index * 3 + 2]*0.05;
                }

                jacob = computeJacob(numvar, numcon, MS, nodes, x);
                var res = computeResidual(numvar, numcon, MS, nodes, x);
                //System.Windows.Forms.MessageBox.Show(res.Norm().ToString());
                S = jacob.Multiply(jacob.T);
                //System.Windows.Forms.MessageBox.Show(S.Det().ToString());
                solve = new ShoNS.Array.SVD(S);
                //if (solve.Rank() < numcon) break;
                //System.Windows.MessageBox.Show(solve.Rank().ToString());
                var dx = jacob.T.Multiply(solve.Solve(-res));
                lN = new List<Line>();
                listNormal.Add(lN);
                foreach (var v in MS.vertices)
                {
                    int index = MS.vertices.IndexOf(v);
                    x[index * 3 + 0, 0] += dx[index * 3 + 0, 0];
                    x[index * 3 + 1, 0] += dx[index * 3 + 1, 0];
                    x[index * 3 + 2, 0] += dx[index * 3 + 2, 0];
                    Vector3d N = new Vector3d();
                    N.X = x[index * 3 + 0, 0];
                    N.Y = x[index * 3 + 1, 0];
                    N.Z = x[index * 3 + 2, 0];
                    lN.Add(new Line(nodes[v.N], nodes[v.N] + N));
                }
            }
        }
    }
}
