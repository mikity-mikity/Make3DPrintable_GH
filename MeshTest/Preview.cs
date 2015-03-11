using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mikity.ghComponents
{
    public partial class MeshTest : Grasshopper.Kernel.GH_Component
    {
        public Rhino.Geometry.Transform zDown_eq = Rhino.Geometry.Transform.Translation(0, 0, 25d);
        public Rhino.Geometry.Transform zDown = Rhino.Geometry.Transform.Translation(0, 0, 15d);
        public Rhino.Geometry.Transform zScale = Rhino.Geometry.Transform.Scale(Rhino.Geometry.Plane.WorldXY, 1, 1, 1d);
        public override void DrawViewportWires(Grasshopper.Kernel.IGH_PreviewArgs args)
        {
            if (Hidden)
            {
                return;
            }
            if (inMesh != null)
            {
                args.Display.DrawMeshWires(newMesh, System.Drawing.Color.Brown);
            }
/*            if (inMesh != null)
            {
                args.Display.DrawMeshShaded(inMesh, new Rhino.Display.DisplayMaterial(System.Drawing.Color.Red));
            }
            if (outMesh != null)
            {
                args.Display.DrawMeshShaded(outMesh, new Rhino.Display.DisplayMaterial(System.Drawing.Color.Red));
            }*/
            if (listNormal != null)
            {
                args.Display.DrawLines(listNormal, System.Drawing.Color.Green);
            }
        }
    }
}
