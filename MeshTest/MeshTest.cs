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
        public double thickness = 1d;
        public List<Line> listNormal;
        private void computeGradient(DoubleArray grad,vertex v,MeshStructure MS, DoubleArray x, Vector3d[] normalPerFaces)
        {
            grad[0, 0] = 0;
            grad[1, 0] = 0;
            grad[2, 0] = 0;
            foreach (var e in v.onering)
            {
                var N = normalPerFaces[e.owner.N];
                int index = MS.vertices.IndexOf(v);
                var cx = x[index * 3 + 0, 0];
                var cy = x[index * 3 + 1, 0];
                var cz = x[index * 3 + 2, 0];
                grad[0,0] += 2 * N.X * N.X * cx + 2 * N.X * N.Y * cy + 2 * N.X * N.Z * cz - 2 * N.X;
                grad[1, 0] += 2 * N.Y * N.Y * cy + 2 * N.Y * N.Z * cz + 2 * N.Y * N.X * cx - 2 * N.Y;
                grad[2, 0] += 2 * N.Z * N.Z * cz + 2 * N.Z * N.X * cx + 2 * N.Z * N.Y * cy - 2 * N.Z;
                var dot = cx * N.X + cy * N.Y + cz * N.Z;
                var norm = N.SquareLength;
                grad[0, 0] += 0.01 * (2 * cx - 4 * dot * N.X + norm * 2 * dot * N.X);
                grad[1, 0] += 0.01 * (2 * cy - 4 * dot * N.Y + norm * 2 * dot * N.Y);
                grad[2, 0] += 0.01 * (2 * cz - 4 * dot * N.Z + norm * 2 * dot * N.Z);
            }
        }
        private void computeHessian(DoubleArray hess, vertex v, MeshStructure MS, DoubleArray x, Vector3d[] normalPerFaces)
        {
            hess[0, 0] = 0;
            hess[0, 1] = 0;
            hess[0, 2] = 0;
            hess[1, 0] = 0;
            hess[1, 1] = 0;
            hess[1, 2] = 0;
            hess[2, 0] = 0;
            hess[2, 1] = 0;
            hess[2, 2] = 0;
            foreach (var e in v.onering)
            {
                var N = normalPerFaces[e.owner.N];
                int index = MS.vertices.IndexOf(v);
                var cx = x[index * 3 + 0, 0];
                var cy = x[index * 3 + 1, 0];
                var cz = x[index * 3 + 2, 0];
                hess[ 0,0] += 2 * N.X * N.X;
                hess[ 1,0] += 2 * N.X * N.Y;
                hess[ 2, 0] += 2 * N.X * N.Z;
                hess[0, 1] += 2 * N.Y * N.X;
                hess[ 1, 1] += 2 * N.Y * N.Y;
                hess[ 2,1] += 2 * N.Y * N.Z;
                hess[0,  2] += 2 * N.Z * N.X;
                hess[ 1,2] += 2 * N.Z * N.Y;
                hess[2, 2] += 2 * N.Z * N.Z;
                var dot = cx * N.X + cy * N.Y + cz * N.Z;
                var norm = N.SquareLength;
                hess[ 0,  0] +=0.01* (2 - 4 * N.X * N.X + norm * 2 * N.X * N.X);
                hess[1, 0] += 0.01 * (-4 * N.Y * N.X + norm * 2 * N.Y * N.X);
                hess[2, 0] += 0.01 * (-4 * N.Z * N.X + norm * 2 * N.Z * N.X);
                hess[1, 1] += 0.01 * (2 - 4 * N.Y * N.Y + norm * 2 * N.Y * N.Y);
                hess[2, 1] += 0.01 * (-4 * N.Z * N.Y + norm * 2 * N.Z * N.Y);
                hess[0, 1] += 0.01 * (-4 * N.X * N.Y + norm * 2 * N.X * N.Y);
                hess[2, 2] += 0.01 * (2 - 4 * N.Z * N.Z + norm * 2 * N.Z * N.Z);
                hess[0, 2] += 0.01 * (-4 * N.X * N.Z + norm * 2 * N.X * N.Z);
                hess[1, 2] += 0.01 * (-4 * N.Y * N.Z + norm * 2 * N.Y * N.Z);
            }
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
            pManager.AddNumberParameter("t", "t", "t", Grasshopper.Kernel.GH_ParamAccess.item);
        }
        protected override void RegisterOutputParams(Grasshopper.Kernel.GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("grad", "grad", "grad", Grasshopper.Kernel.GH_ParamAccess.list);
        }
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            Mesh myMesh=new Mesh();
            if (!DA.GetData(0, ref myMesh)) { return; }
            if (!DA.GetData(1, ref thickness)) { return; }
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
            DoubleArray q = DoubleArray.Zeros(numvar, 1);
            List<Point3d> nodes = new List<Point3d>();
            Vector3d[] normalPerFaces = new Vector3d[MS.nFaces];
            nodes.AddRange(myMesh.Vertices.ToPoint3dArray());

            //initial guess
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
                N.Unitize();
                x[index * 3 + 0, 0] = N[0];
                x[index * 3 + 1, 0] = N[1];
                x[index * 3 + 2, 0] = N[2];
            }
            List<double> grads = new List<double>();
            List<double> residuals = new List<double>();
            double tol = 0.00000001;
            var grad = new DoubleArray(3, 1);
            var hess = new DoubleArray(3, 3);
            foreach (var v in MS.vertices)
            {
                for (int i = 0; i < 500; i++)
                {
                    computeGradient(grad, v, MS, x, normalPerFaces);
                    computeHessian(hess, v, MS, x, normalPerFaces);
                    var sol = new LU(hess);
                    var dx = sol.Solve(-grad);
                    int index = MS.vertices.IndexOf(v);
                    x[index * 3 + 0, 0] += dx[0, 0];
                    x[index * 3 + 1, 0] += dx[1, 0];
                    x[index * 3 + 2, 0] += dx[2, 0];
                    if (grad.Norm() < tol)
                    {
                        grads.Add(i);
                        break;
                    }
                }
            }
            inMesh = myMesh.DuplicateMesh();
            outMesh = myMesh.DuplicateMesh();
            listNormal = new List<Line>();
            foreach (var v in MS.vertices)
            {
                int index=MS.vertices.IndexOf(v);
                inMesh.Vertices[v.N] = new Point3f((float)(nodes[v.N].X + x[index * 3 + 0, 0] * thickness / 2d), (float)(nodes[v.N].Y + x[index * 3 + 1, 0] * thickness / 2d), (float)(nodes[v.N].Z + x[index * 3 + 2, 0] * thickness / 2d));
                outMesh.Vertices[v.N] = new Point3f((float)(nodes[v.N].X - x[index * 3 + 0, 0] * thickness / 2d), (float)(nodes[v.N].Y - x[index * 3 + 1, 0] * thickness / 2d), (float)(nodes[v.N].Z - x[index * 3 + 2, 0] * thickness / 2d));
                Vector3d N = new Vector3d();
                N.X = x[index * 3 + 0, 0];
                N.Y = x[index * 3 + 1, 0];
                N.Z = x[index * 3 + 2, 0];

                listNormal.Add(new Line(nodes[v.N], nodes[v.N] + N));
            }
            DA.SetDataList(0, grads);
        }
    }
}
