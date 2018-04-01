using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;

namespace GalaxyGen
{
    public class Consts
    {
        //Begining of user defined inputs
        /// <summary>
        /// self-explanitory
        /// </summary>
        public const double Radius = 256;

        public const int ResX = 588;
        public const int ResY = 624;

        /// <summary>
        /// The resoultion of the grid
        /// </summary>
		public const int Resolution = 16;

        /// <summary>
        /// horizontal subdivisions inside each gridbox
        /// </summary>
		public const int HsubDivisions = 64;

        /// <summary>
        /// vertical subdivisions inside each gridbox
        /// </summary>
		public const int VsubDivisions = 16;

        //End of user defined inputs

        public const double Tau = Math.PI * 2;

        public const double VDivision = (Radius / Resolution);

        public const double VsDivision = (VDivision / VsubDivisions);

        public const double HDivision = (Tau / Resolution);

        public const double HsDivision = (HDivision / HsubDivisions);

    }

    public class Rendering
    {
        private static PaintEventArgs e;

        private static int[] defaultcolour = { 255, 0, 0 };

        public static void Init(PaintEventArgs E)
        {
            e = E;
        }

        public static void Flip()
        {
            SolidBrush b = new SolidBrush(Color.Black);
            e.Graphics.FillRectangle(b, 0, 0, Consts.ResX, Consts.ResY);
        }
    }

    //Perlin noise class from https://lotsacode.wordpress.com/2010/02/24/perlin-noise-in-c/
    public class PerlinNoise
    {
        private const int GradientSizeTable = 256;
        private readonly Random _random;
        private readonly double[] _gradients = new double[GradientSizeTable * 3];
        /* Borrowed from Darwyn Peachey (see references above).
           The gradient table is indexed with an XYZ triplet, which is first turned
           into a single random index using a lookup in this table. The table simply
           contains all numbers in [0..255] in random order. */
        private readonly byte[] _perm = new byte[] {
              225,155,210,108,175,199,221,144,203,116, 70,213, 69,158, 33,252,
                5, 82,173,133,222,139,174, 27,  9, 71, 90,246, 75,130, 91,191,
              169,138,  2,151,194,235, 81,  7, 25,113,228,159,205,253,134,142,
              248, 65,224,217, 22,121,229, 63, 89,103, 96,104,156, 17,201,129,
               36,  8,165,110,237,117,231, 56,132,211,152, 20,181,111,239,218,
              170,163, 51,172,157, 47, 80,212,176,250, 87, 49, 99,242,136,189,
              162,115, 44, 43,124, 94,150, 16,141,247, 32, 10,198,223,255, 72,
               53,131, 84, 57,220,197, 58, 50,208, 11,241, 28,  3,192, 62,202,
               18,215,153, 24, 76, 41, 15,179, 39, 46, 55,  6,128,167, 23,188,
              106, 34,187,140,164, 73,112,182,244,195,227, 13, 35, 77,196,185,
               26,200,226,119, 31,123,168,125,249, 68,183,230,177,135,160,180,
               12,  1,243,148,102,166, 38,238,251, 37,240,126, 64, 74,161, 40,
              184,149,171,178,101, 66, 29, 59,146, 61,254,107, 42, 86,154,  4,
              236,232,120, 21,233,209, 45, 98,193,114, 78, 19,206, 14,118,127,
               48, 79,147, 85, 30,207,219, 54, 88,234,190,122, 95, 67,143,109,
              137,214,145, 93, 92,100,245,  0,216,186, 60, 83,105, 97,204, 52};

        public PerlinNoise(Random r)
        {
            _random = r;
            InitGradients();
        }

        public double Noise(double x, double y, double z)
        {
            /* The main noise function. Looks up the pseudorandom gradients at the nearest
               lattice points, dots them with the input vector, and interpolates the
               results to produce a single output value in [0, 1] range. */

            int ix = (int)Math.Floor(x);
            double fx0 = x - ix;
            double fx1 = fx0 - 1;
            double wx = Smooth(fx0);

            int iy = (int)Math.Floor(y);
            double fy0 = y - iy;
            double fy1 = fy0 - 1;
            double wy = Smooth(fy0);

            int iz = (int)Math.Floor(z);
            double fz0 = z - iz;
            double fz1 = fz0 - 1;
            double wz = Smooth(fz0);

            double vx0 = Lattice(ix, iy, iz, fx0, fy0, fz0);
            double vx1 = Lattice(ix + 1, iy, iz, fx1, fy0, fz0);
            double vy0 = Lerp(wx, vx0, vx1);

            vx0 = Lattice(ix, iy + 1, iz, fx0, fy1, fz0);
            vx1 = Lattice(ix + 1, iy + 1, iz, fx1, fy1, fz0);
            double vy1 = Lerp(wx, vx0, vx1);

            double vz0 = Lerp(wy, vy0, vy1);

            vx0 = Lattice(ix, iy, iz + 1, fx0, fy0, fz1);
            vx1 = Lattice(ix + 1, iy, iz + 1, fx1, fy0, fz1);
            vy0 = Lerp(wx, vx0, vx1);

            vx0 = Lattice(ix, iy + 1, iz + 1, fx0, fy1, fz1);
            vx1 = Lattice(ix + 1, iy + 1, iz + 1, fx1, fy1, fz1);
            vy1 = Lerp(wx, vx0, vx1);

            double vz1 = Lerp(wy, vy0, vy1);
            return Lerp(wz, vz0, vz1);
        }

        private void InitGradients()
        {
            for (int i = 0; i < GradientSizeTable; i++)
            {
                double z = 1f - 2f * _random.NextDouble();
                double r = Math.Sqrt(1f - z * z);
                double theta = 2 * Math.PI * _random.NextDouble();
                _gradients[i * 3] = r * Math.Cos(theta);
                _gradients[i * 3 + 1] = r * Math.Sin(theta);
                _gradients[i * 3 + 2] = z;
            }
        }

        private int Permutate(int x)
        {
            const int mask = GradientSizeTable - 1;
            return _perm[x & mask];
        }

        private int Index(int ix, int iy, int iz)
        {
            // Turn an XYZ triplet into a single gradient table index.
            return Permutate(ix + Permutate(iy + Permutate(iz)));
        }

        private double Lattice(int ix, int iy, int iz, double fx, double fy, double fz)
        {
            // Look up a random gradient at [ix,iy,iz] and dot it with the [fx,fy,fz] vector.
            int index = Index(ix, iy, iz);
            int g = index * 3;
            return _gradients[g] * fx + _gradients[g + 1] * fy + _gradients[g + 2] * fz;
        }

        private double Lerp(double t, double value0, double value1)
        {
            // Simple linear interpolation.
            return value0 + t * (value1 - value0);
        }

        private double Smooth(double x)
        {
            /* Smoothing curve. This is used to calculate interpolants so that the noise
              doesn't look blocky when the frequency is low. */
            return x * x * (3 - 2 * x);
        }
    }

    public class Debug
    {
        /// <summary>
        /// boolean to enable/disable the debug functions
        /// </summary>
        public const bool Active = true;

        //stack used to hold the timestamps of every check, so that multiple one can be queued at once
        private static Stack<int> timestamps = new Stack<int>();

        //needed to paint stuff the the screen
        private static PaintEventArgs e;

        //the default colour used by Gplot if none is passed as an argument
        private static int[] defaultcolour = { 255, 0, 0 };

        /// <summary>
        /// Collect the event args needed for drawing to the canvas
        /// </summary>
        /// <param name="E">PaintEventArgs needed for drawing</param>
        public static void InitG(PaintEventArgs E)
        {
            e = E;
        }


        /// <summary>
        /// Draw a dot to the canvas
        /// </summary>
        /// <param name="point">the position to draw to</param>
        /// <param name="size">the size of the dot</param>
        /// <param name="rgb">the colour of the dot</param>
        public static void Gplot(RealPoint point, double size = 1, int[] rgb = null)
        {
            rgb = rgb ?? defaultcolour; //if rgb is null, use the default colour
            SolidBrush p = new System.Drawing.SolidBrush(Color.FromArgb(rgb[0], rgb[1], rgb[2])); //create a brush and convert the colour from int[]

            //Draw a circle
            e.Graphics.FillEllipse(p, new RectangleF(new PointF((float)(point.x - (size / 2)), (float)(point.y - (size / 2))), new SizeF((float)size, (float)size)));
        }

        public static void Gwrite(RealPoint point, string msg, int[] rgb = null, params object[] args)
        {
            rgb = rgb ?? defaultcolour; //if rgb is null, use the default colour
            SolidBrush p = new System.Drawing.SolidBrush(Color.FromArgb(rgb[0], rgb[1], rgb[2])); //create a brush and convert the colour from int[]
            
            //Draw a circle
            e.Graphics.DrawString(string.Format(msg, args),new Font("Arial", 16),p,(float)point.x, (float)point.y);
        }

        /// <summary>
        /// queue up a timestamp
        /// </summary>
        public static void check()
        {
            if (Active)
            {
                timestamps.Push(Environment.TickCount);
            }
        }

        /// <summary>
        /// A combinination of both end() and check()
        /// </summary>
        /// <returns>the amount of time elapsed since the timestamp at the top of the stack, in ticks</returns>
        public static string reset()
        {
            if (Active)
            {
                int temp = timestamps.Pop();
                timestamps.Push(Environment.TickCount);
                return (Environment.TickCount - temp) + " ticks";
            }
        }

        /// <summary>
        /// pulls the last timestamp from the stack
        /// </summary>
        /// <returns>the amount of time elapsed since the timestamp at the top of the stack, in ticks</returns>
        public static string end()
        {
            if (Active)
            {
                return (Environment.TickCount - timestamps.Pop()) + " ticks";
            }
        }

        /// <summary>
        /// Print
        /// </summary>
        public static void p(string Word, params object[] args) { Console.WriteLine(Word, args); }
    }

    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            this.Width = Consts.ResX;
            this.Height = Consts.ResY;
            this.BackColor = Color.Black;
        }

        private void Form1_Paint(object sender, PaintEventArgs e)
        {
            //Main
            Debug.InitG(e);
            Rendering.Init(e);
            Galaxy g = new Galaxy(1616568269);
        }
    }

    //Generic abstract classes begin

    /// <summary>
    /// Class used to repesent a point in 2d space
    /// </summary>
    public class RealPoint
    {
        public double x;
        public double y;

        public RealPoint(double X, double Y)
        {
            x = X;
            y = Y;
        }

        //Secondary constructor used as a substitute for the now-defunct referencepoint class
        public RealPoint(double Theta, double R, RealPoint Centre)
        {
            x = Centre.x + (R * Math.Cos(Theta));
            y = Centre.y + (R * Math.Sin(Theta));
        }

        //extra constuctor to allow creation of empty points
        public RealPoint()
        {
            x = 0;
            y = 0;
        }

        /// <summary>
        /// Get a point on a line between this point and another
        /// </summary>
        /// <param name="p2">the second point</param>
        /// <param name="a">the position on the line, between 1 and 0. 1 == p2 </param>
        public RealPoint Lerp(RealPoint p2, double a)
        {
            if (a < 1 && a > 0)
            {
                return (this * (1 - a)) + (p2 * a);
            }
            else
            {
                return (a >= 1) ? p2 : this;
            }
        }

        /// <summary>
        /// Move this point along a line between the current position and another point
        /// </summary>
        /// <param name="p2">the second point</param>
        /// <param name="a">the position on the line, between 1 and 0. 1 == p2 </param>
        public void SetLerp(RealPoint p2, double a)
        {
            RealPoint b = Lerp(p2, a);
            x = b.x; y = b.y;
        }

        /// <summary>
        /// Get the gradient of a line drawn between this point and another
        /// </summary>
        /// <param name="p2">the second point</param>
        public double Grad(RealPoint p2)
        {
            return (p2.y - y) / (p2.x - x);
        }

        /// <summary>
        /// Get the distance between two points
        /// </summary>
        /// <param name="p1">the first point</param>
        /// <param name="p2">the second point</param>
        /// <returns></returns>
        public static double Distance(RealPoint p1, RealPoint p2)
        {
            return Math.Sqrt((Math.Abs(p1.x - p2.x) * Math.Abs(p1.x - p2.x)) + Math.Abs(p1.y - p2.y) * Math.Abs(p1.y - p2.y));
        }

        public static RealPoint operator +(RealPoint p1, RealPoint p2)
        {
            return new RealPoint(p1.x + p2.x, p1.y + p2.y);
        }
        public static RealPoint operator -(RealPoint p1, RealPoint p2)
        {
            return new RealPoint(p1.x - p2.x, p1.y - p2.y);
        }
        public static RealPoint operator *(RealPoint p1, RealPoint p2)
        {
            return new RealPoint(p1.x * p2.x, p1.y * p2.y);
        }
        public static RealPoint operator /(RealPoint p1, RealPoint p2)
        {
            return new RealPoint(p1.x / p2.x, p1.y / p2.y);
        }
        public static RealPoint operator *(RealPoint p1, double a)
        {
            return new RealPoint(p1.x * a, p1.y * a);
        }
        public static RealPoint operator /(RealPoint p1, double a)
        {
            return new RealPoint(p1.x / a, p1.y / a);
        }
        public static RealPoint operator *(double a, RealPoint p1)
        {
            return new RealPoint(p1.x * a, p1.y * a);
        }
        public static RealPoint operator /(double a, RealPoint p1)
        {
            return new RealPoint(p1.x / a, p1.y / a);
        }
        public override string ToString()
        {
            return x + ", " + y;
        }
        public bool Equals(RealPoint p2)
        {
            return x == p2.x && p2.y == y;
        }
    }

    public class Line
    {
        // Y = Xmod(x) + Ymod
        public double Xmod;
        public double Ymod;

        //the two points used to create the line
        public RealPoint p1;
        public RealPoint p2;

        public Line(RealPoint start, RealPoint end)
        {
            Xmod = start.Grad(end);
            Ymod = -((Xmod * start.x)-start.y);

            p1 = start;
            p2 = end;
        }

        public Line(RealPoint point, double Grad)
        {
            Xmod = Grad;
            Ymod = -((Xmod * point.x) - point.y);

            p1 = point;
            p2 = point;
        }

        public RealPoint Intersect(Line other)
        {
            if (other.Xmod == Xmod && other.Ymod == Ymod)
            {
                return null;
            }

            double x = (other.Ymod - Ymod) / (Xmod - other.Xmod);

            return new RealPoint(x, Xmod * (x) + Ymod);
        }

        public Line ShiftPerp(double shift)
        {
            Debug.p("creating shift with distance of {0} from target of {1}", RealPoint.Distance(new RealPoint(Math.Atan(Xmod) + (0.5 * Math.PI), shift, p1), p1), shift);
            return new Line(new RealPoint(Math.Atan(Xmod)+(0.5*Math.PI), shift,p1), Xmod);
        }
        
        public double f(double x)
        {
            return (Xmod * x) + Ymod;
        }

        public void Draw(int[] Colour = null, int start = 0, int end = Consts.ResX)
        {
            for (int i = start; i < end; i++)
            {
                Debug.Gplot(new RealPoint(i,f(i)), 1, Colour);
            }
        }
    }

    public abstract class Default
    {
        public dynamic Parent;

        public Default(dynamic parent)
        {
            Parent = parent;
        }
    }


    public abstract class Arc : Default
    {
        /// <summary>
        /// The height of the top radius, in VsubDivisions
        /// </summary>
        public int TopR;
        /// <summary>
        /// The height of the top radius, in VsubDivisions
        /// </summary>
        public int BottomR;
        /// <summary>
        /// The theta of the start of the arc, in HsubDivision
        /// </summary>
        public int StartTheta;
        /// <summary>
        /// The theta of the end of the arc, in HsubDivisions
        /// </summary>
        public int EndTheta;

        public Arc(int aR, int bR, int aT, int bT, Default parent) : base(parent)
        {
            TopR = aR;
            BottomR = bR;
            StartTheta = aT;
            EndTheta = bT;
        }

        public void Draw(RealPoint Centre)
        {
            for (double i = 0; i < 1; i += 0.01)
            {
                Debug.Gplot(new RealPoint((StartTheta * Consts.HsDivision) + (((EndTheta * Consts.HsDivision) - (StartTheta * Consts.HsDivision)) * i), (TopR * Consts.VsDivision), Centre), 1, new int[3] { 0, 0, 255 });
                Debug.Gplot(new RealPoint((StartTheta * Consts.HsDivision) + (((EndTheta * Consts.HsDivision) - (StartTheta * Consts.HsDivision)) * i), (BottomR * Consts.VsDivision), Centre), 1, new int[3] { 0, 0, 255 });
                Debug.Gplot(new RealPoint((StartTheta * Consts.HsDivision), (BottomR * Consts.VsDivision) + (((TopR * Consts.VsDivision) - (BottomR * Consts.VsDivision)) * i), Centre));
                Debug.Gplot(new RealPoint((EndTheta * Consts.HsDivision), (BottomR * Consts.VsDivision) + (((TopR * Consts.VsDivision) - (BottomR * Consts.VsDivision)) * i), Centre));
            }
        }
    }

    //end abstract classes

    //begin base classes

    public class Galaxy : Default
    {
        public int Seed;
        public Random r { get; private set; }

        public Galaxy(int seed) : base(null)
        {
            Seed = seed;

            Main();
        }

        public Galaxy(string seed) : base(null)
        {
            string TrueSeed = "";
            foreach (char i in seed)
            {
                TrueSeed += ((int)i * 69).ToString();
            }
            Console.WriteLine(TrueSeed);
            Seed = TrueSeed.GetHashCode();

            Main();
        }

        public Galaxy() : base(null)
        {
            Seed = new Random().Next(0, Int32.MaxValue);

            Main();
        }

        private void Main()
        {
            for (double i = 0; i < Math.PI * 2; i += 0.001)
            {
                Debug.Gplot(new RealPoint(i, Consts.Radius, Centre));
            }
            for (double i = 0; i < Math.PI * 2; i += 0.01)
            {
                Debug.Gplot(new RealPoint(i, CentreGap * Consts.VDivision, Centre));
            }

            PerlinNoise n = new PerlinNoise(new Random());

            GenHLs();
            GenSectors(32,64,20,0.4,1,1,10);
            //TODO: Add function for bites and plan system for planet distribution
        }

        private void GenHLs()
        {
            r = new Random(Seed);
            HLs = new Hyperlane[4];
            for (int i = 0; i < 4; i++)
            {
                HLs[i] = new Hyperlane(3, 200, 0.002, this);
            }
        }

        public void DragToClosestHL(SolarSystem ToDrag)
        {
            RealPoint Closest = ClosestHLPoint(ToDrag.Position);
            double dist = RealPoint.Distance(ToDrag.Position, Closest);

            //move the point a percentage distance towards the lane specified by the curve -0.01x^2 + 1
            ToDrag.Position.SetLerp(Closest, Math.Max((-0.01 * (dist * dist)) + 1, 0));
        }

        public RealPoint ClosestHLPoint(RealPoint From)
        {
            double lowest = double.MaxValue; //Set the lowest distance found to a distance that's impossible to beat
            RealPoint Closest = new RealPoint(); //Generate a new point to stop compiler errors
            foreach (Hyperlane i in HLs) //check through each hyperlane
            {
                //get the closest point on the current hyperlane
                RealPoint checkpoint = i.GetAt(From);
                double checkdistance = RealPoint.Distance(checkpoint, From);

                //check if it is closer than the previous closest point
                if (checkdistance < lowest)
                {
                    lowest = checkdistance; //use the previously generated point and distance for efficency
                    Closest = checkpoint;
                }
            }
            return Closest;
        }

        public double ClosestHLDist(RealPoint From)
        {
            double lowest = double.MaxValue; //Set the lowest distance found to a distance that's impossible to beat
            RealPoint Closest = new RealPoint(); //Generate a new point to stop compiler errors
            foreach (Hyperlane i in HLs) //check through each hyperlane
            {
                //get the closest point on the current hyperlane
                RealPoint checkpoint = i.GetAt(From);
                double checkdistance = RealPoint.Distance(checkpoint, From);

                //check if it is closer than the previous closest point
                if (checkdistance < lowest)
                {
                    lowest = checkdistance; //use the previously generated point and distance for efficency
                    Closest = checkpoint;
                }
            }
            return lowest;
        }

        /// <param name="MinSize">Minimum size of a sector, in HsubDivisions</param>
        /// <param name="MaxSize">Maximum size of a sector, in HsubDivisions</param>
        /// <param name="MinHeight">Minimum height of a ???</param>
        /// <param name="MaxHeight"></param>
        private void GenSectors(int MinLength, int MaxLength, int BiteProb, double MinBiteLength, double MaxBiteLength, int MinBiteHeight, int MaxBiteHeight)
        {
            r = new Random(Seed);
            Sectors = new Run[Consts.Resolution - CentreGap];
            for (int i = Consts.Resolution-CentreGap; i > 0; i--)
            {
                Sectors[i-1] = new Run(i+CentreGap, (Consts.HsubDivisions*Consts.Resolution)/MaxLength, (Consts.HsubDivisions * Consts.Resolution) / MinLength, 0.5, this);
                foreach (Sector s in Sectors[i - 1].Sectors)
                {
                    if (r.Next(0, 100) <= BiteProb && !s.CheckHLCollision(HLs))
                    {
                        s.Break(MinBiteLength, MaxBiteLength, MinBiteHeight, MaxBiteHeight);
                    }
                }
            }
        }

        Hyperlane[] HLs;
        public int CentreGap = 2;

        Run[] Sectors;

        public RealPoint Centre = new RealPoint(Consts.Radius, Consts.Radius);
    }

    public class Hyperlane : Default
    {
        //The line to be offset by perlin noise
        public Line Base;
        
        PerlinNoise noise; //an instance of the perlin noise class, used to get the noise
        double frequency = 0.5; //the frequency of the noise
        double Y; //Since perlin noise requires an X, Y and Z, the Y to be used is stored here

        public double MaxSize; //the amount to muliply the perlin noise by, and hence the largest offset
        
        /// <param name="MaxLineDiversion">The maximum offset for the line from going directly through the centre</param>
        /// <param name="MaxExpanse">the maximum amount the perlin noise can deviate the base line</param>
        /// <param name="Frequency">the frequency of the perlin noise</param>
        public Hyperlane(int MaxLineDiversion, double MaxExpanse, double Frequency, Object parent) : base(parent)
        {
            Random R = Parent.r; //Get the parent's random function, if we generated it here every lane would be the same

            Y = R.Next(); //Generate the Y component
            MaxSize = MaxExpanse; //Set the max offset
            frequency = Frequency; //Set the frequency

            //Get the theta of the first point of the line by picking one of the HsubDivisions
            double InitialTheta = R.Next(0, Consts.Resolution * Consts.HsubDivisions) * Consts.HsDivision;
            
            double Diversion = 
                (Parent.CentreGap * Consts.VDivision) +  //Get the real size of the centre by converting it to real units
                Math.Max(R.Next(0,MaxLineDiversion) * Consts.HsDivision,  //Add it to either a random number between 0 and the maximum offset
                (MaxExpanse*0.3)); //or the amplitude of the perlin noise, so that it's less likely to extend into the centre

            //The two points of the chord are created by projecting out from a point around the centre
            RealPoint Pivot = new RealPoint(InitialTheta, Diversion, Parent.Centre);

            //Calculate the length using pythag (it's correct I've checked)
            double ChordLen = Math.Sqrt((-Diversion * -Diversion) + (Consts.Radius * Consts.Radius));

            //Generate the bisecting line of the chord
            Base = new Line( new RealPoint(InitialTheta+(Math.PI*0.5), ChordLen, Pivot), new RealPoint(InitialTheta - (Math.PI * 0.5), ChordLen, Pivot));

            //Generate a new instance of the perlin noise class
            noise = new PerlinNoise(R);

            Base.Draw();
            Base.ShiftPerp(MaxSize).Draw();
            Base.ShiftPerp(-MaxSize).Draw();
            for (double i = 0; i < 1; i += 0.005)
            {
                Debug.Gplot(GetAt(Base.p1.Lerp(Base.p2,i)));
            }
            Rendering.Flip();
        }

        /// <summary>
        /// Get a close-ish point on the hyperlane to a given point 
        /// </summary>
        public RealPoint GetAt(RealPoint Point)
        {
            //Get the intercept point between the base line and another line at a right angle that passes through the point
            RealPoint i = Base.Intersect(new Line(Point,-Base.Xmod));

            //Get the perlin noise at the location, modified by the frequency so that it is correct
            double n = noise.Noise(i.x*frequency, Y*frequency, i.y*frequency);

            return new RealPoint(
                Math.Atan2(Base.p2.y - Base.p1.y, Base.p2.x - Base.p1.x) //convert the line's gradient to radians
                + (0.5d * Math.PI), //turn it 90 degrees so it's at a right angle
                n * MaxSize, //project out by the noise
                i //from the intercept point
            );
        }

        public double GetAtRAW(RealPoint Point)
        {
            //Get the intercept point between the base line and another line at a right angle that passes through the point
            RealPoint i = Base.Intersect(new Line(Point, -Base.Xmod));
            
            return noise.Noise(i.x * frequency, Y * frequency, i.y * frequency);
        }
    }

    public class Run : Default
    {
        int Radius;
        public Sector[] Sectors;

        public Run(int r, int MinSectors, int MaxSectors, double core, Galaxy parent) : base(parent)
        {
            int SectorCount = parent.r.Next(MinSectors, MaxSectors);
            Sectors = new Sector[(SectorCount % 2 == 1) ? SectorCount-1 : SectorCount];

            int shift = parent.r.Next(0,Consts.Resolution*Consts.HsubDivisions);
            int Base = (int)Math.Floor((double)((Consts.Resolution * Consts.HsubDivisions) / Sectors.Length));
            int Remainder = (Base - (int)Math.Ceiling(Base * core));

            int[] lengths = new int[Sectors.Length];
            for (int i = 0; i < Sectors.Length; i++)
            {
                //loop to precalculate all the length, so we don't waste resources regenerating the first sector
                //takes more memory due to having everything at once rather than overwriting one variable repeatedly
                //but that's a sacrifice I'm willing to make
                lengths[i] = parent.r.Next(0,Remainder);
            }
            
            //Set the first sector seperately because you can't do index [-1] for the last element like in python
            Sectors[0] = new Sector(
                                            r * Consts.VsubDivisions, //set the radius, self explanitory
                                            (r - 1) * Consts.VsubDivisions, //set the radius, self explanitory
                                            shift - lengths[lengths.Length-1] - (Consts.Resolution * Consts.HsubDivisions) % Sectors.Length, //take up the space left by the last guy
                                            shift + Base - lengths[0], //expand by the base amount, minus the next guy's extention
                                            this //parent
                                            );

            for (int i = 1; i < Sectors.Length; i++)
            {
                shift += Base;
                Sectors[i] = new Sector(
                                            r*Consts.VsubDivisions, //set the radius, self explanitory
                                            (r-1)*Consts.VsubDivisions, //set the radius, self explanitory
                                            shift - lengths[i-1], //take up the space left by the last guy
                                            shift + Base - lengths[i], //expand by the base amount, minus the next guy's extention
                                            this //parent
                                            );
            }
        }
    }

    public class Sector : Arc
    {
        public Sector(int TopRadius, int BottomRadius, int Start, int End, Default Parent) : base(TopRadius, BottomRadius, Start, End, Parent)
        {
            Subsectors = new Subsector[(TopRadius-BottomRadius), (End - Start)];
            int yindex = 0;
            int xindex = 0;
            for (int y = TopRadius; y > BottomRadius; y--)
            {
                for (int x = StartTheta; x < EndTheta; x++)
                {
                    Subsectors[yindex, xindex] = new Subsector(y, y-1, x, x+1, this);
                    xindex++;
                }
                yindex++;
                xindex = 0;
            }
            //drawing
            //Draw(Parent.Parent.Centre);
        }

        public void Break(double MinLength, double MaxLength, int MinHeight, int MaxHeight)
        {
            int Length = Parent.Parent.r.Next((int)(MinLength*Subsectors.GetLength(1)), (int)(MaxLength*Subsectors.GetLength(1)));
            int From = Parent.Parent.r.Next(0, Subsectors.GetLength(1) - 1 - Length);
            int Height = Parent.Parent.r.Next(MinHeight, MaxHeight);
            int Extention = (Parent.Parent.r.Next(0, 2) == 0) ? 0 : (Subsectors.GetLength(0) - 1 - Height);
            for (int x = From; x < From+Length; x++)
            {
                for (int y = 0; y < Height; y++)
                {
                    Subsectors[y + Extention, x].SysMult = 0;
                    //Subsectors[y + Extention, x].Draw(Parent.Parent.Centre);
                }
            }
            
            new Subsector(Subsectors[Extention, From].TopR, Subsectors[Height+Extention, From+Length].BottomR, Subsectors[Extention, From].StartTheta, Subsectors[Height + Extention, From + Length].EndTheta, this).Draw(Parent.Parent.Centre);
        }

        public bool CheckHLCollision(Hyperlane[] HLlist)
        {
            foreach (Hyperlane i in HLlist)
            {
                if (CheckHLCol(i)) { return true; }
            }
            return false;
        }

        bool CheckHLCol(Hyperlane ColWith)
        {
            Draw(Parent.Parent.Centre);
            
            Line BoundLine1 = new Line(new RealPoint(StartTheta, TopR, Parent.Parent.Centre), new RealPoint(StartTheta, BottomR, Parent.Parent.Centre));
            Line BoundLine2 = new Line(new RealPoint(EndTheta, TopR, Parent.Parent.Centre), new RealPoint(EndTheta, BottomR, Parent.Parent.Centre));
            
            BoundLine1.Draw(new int[3] { 0, 0, 255 });
            BoundLine2.Draw(new int[3] { 0, 0, 255 });
            ColWith.Base.Draw(new int[3] { 0, 255, 0 });
            ColWith.Base.ShiftPerp(ColWith.MaxSize).Draw(new int[3] { 255, 0, 0 });
            ColWith.Base.ShiftPerp(-ColWith.MaxSize).Draw(new int[3] { 255, 0, 0 });
            //Debug.p("Getting distance of {0} with a max size of {1}", RealPoint.Distance(new RealPoint(0,ColWith.Base.ShiftPerp(ColWith.MaxSize).f(0)), new RealPoint(0, ColWith.Base.f(0))), ColWith.MaxSize);
            
            double IntersectR1a = RealPoint.Distance(ColWith.Base.ShiftPerp(ColWith.MaxSize).Intersect(BoundLine1), Parent.Parent.Centre);
            double IntersectR2a = RealPoint.Distance(ColWith.Base.ShiftPerp(-ColWith.MaxSize).Intersect(BoundLine1), Parent.Parent.Centre);
            double IntersectR1b = RealPoint.Distance(ColWith.Base.ShiftPerp(ColWith.MaxSize).Intersect(BoundLine2), Parent.Parent.Centre);
            double IntersectR2b = RealPoint.Distance(ColWith.Base.ShiftPerp(-ColWith.MaxSize).Intersect(BoundLine2), Parent.Parent.Centre);
            if (!(IntersectR1a <= TopR && IntersectR1a >= BottomR)&& !(IntersectR2a <= TopR && IntersectR2a >= BottomR) && !(IntersectR1b <= TopR && IntersectR1b >= BottomR) && !(IntersectR2b <= TopR && IntersectR2b >= BottomR))
            {
                return false;
            }
            return true;
        }

        Subsector[,] Subsectors;
    }

    public class Subsector : Arc
    {
        public Subsector(int TopRadius, int BottomRadius, int Start, int End, Default Parent) : base(TopRadius, BottomRadius, Start, End, Parent) { }

        SolarSystem[] Systems;
        public double SysMult = 1;
    }
    
    public class SolarSystem : Default
    {
        public RealPoint Position;

        SolarSystem(RealPoint a, Default parent) : base(parent)
        {
            Position = a;
            Parent.Parent.Parent.DragToClosestHL(this);
        }
    }

    public abstract class Body : Default
    {
        double Radius;

        public Body(double radius, Default parent) : base(parent) { Radius = radius; }
    }

    public class Sun : Body
    {
        Sun(double radius, Default parent) : base(radius, parent) { }
    }
    
    public class Planet : Body
    {
        object Materials;
        double OrbitalRadius;
        double RotationalPeriod;

        Planet(double radius, double orbit, double rotation, Default parent) : base(radius, parent)
        {
            OrbitalRadius = orbit;
            RotationalPeriod = rotation;
        }
    }

    /*
	public class Cursor
	{

		Pen c;
		Color Default;
		Graphics g;

		public Cursor(PaintEventArgs e, Color d)
		{
			g = e.Graphics;
			Default = d;
			c = new Pen(d);
		}

        public void Dot(float x1, float y1, float Sizex, float Sizey, Color C = new Color())
        {
            Brush b = new SolidBrush(C);
            g.FillEllipse(b, x1, y1, Sizex, Sizey);
        }

		public void Line(float x1, float y1, float x2, float y2, Color C = new Color())
		{
			c = new Pen((C==new Color())?Default:C);
			g.DrawLine(c, x1, y1, x2, y2); 
		}

		public void Line(double x1, double y1, double x2, double y2, Color C = new Color())
		{
			c = new Pen((C == new Color()) ? Default : C);
			g.DrawLine(c, (float)x1, (float)y1, (float)x2, (float)y2);
		}

		public void Arc(int Tradius, int Start, int End, double Vscale = Consts.Division, double precision = 0.05, Color C = new Color())
		{
			for (double i = Start; i < End; i += precision)
			{
				Line(
				Consts.Radius + (Tradius * Vscale) * Math.Cos(i * Consts.Segment),
				Consts.Radius + (Tradius * Vscale) * Math.Sin(i * Consts.Segment),
				Consts.Radius + (Tradius * Vscale) * Math.Cos((i + precision) * Consts.Segment),
				Consts.Radius + (Tradius * Vscale) * Math.Sin((i + precision) * Consts.Segment),
				C
				);
			}
		}
    }

	public class Point
	{
		public double Theta;
		public int r;

		public Point(int R, int T)
		{
			Theta = T;
			r = R;
		}

        public RealPoint ToReal( double CentreX = Consts.Radius, double CentreY = Consts.Radius)
        {
            return new RealPoint( CentreX + (r * Math.Cos(Theta)), CentreY + (r * Math.Sin(Theta)) );
        }
	}
    */


    /*
	public class FilledArc
	{
		public int Tradius;
		public int Start;
		public int End;

		protected FilledArc(int a, int b, int c)
		{
			Tradius = a;
			Start = b;
			End = c;
		}

		protected FilledArc() { }
	}

    public class Spline
    {
        Point Start;
        Point End;

        float[] deviantcies; //lists deviances as a percentage of total line length

        public Spline(Point S, Point E, int D)
        {
            Start = S;
            End = E;
            deviantcies = new float[D];
        }


    }

    public class Galaxy
    {
        Random R;

        public int Hyperlanes = 3;
        public int deviations = 6;
        public RealPoint[,] HLpoints;

        public int MinHLseperation = 8;

        int MaxDeviationSize = 6;

        public Galaxy(Random r)
        {
            HLpoints = new RealPoint[Hyperlanes,deviations+2];
            R = r;
        }

        public void BuildHyperlanes()
        {
            for (int i = 0; i < Hyperlanes; i++)
            {
                int Ot = R.Next(0, Consts.Divisions * Consts.IsubDivisions);
                HLpoints[i, 0] = new Point(Consts.Radius, Ot).ToReal();
                int Tt = Ot + R.Next(MinHLseperation, (Consts.Divisions * Consts.IsubDivisions)*MinHLseperation) % (Consts.Divisions * Consts.IsubDivisions);
                
                HLpoints[i, deviations + 1] = new Point(Consts.Radius, Tt).ToReal();
                for (int p = 0; p < deviations; p++)
                {
                    //nobody toucha my spaget code
                    RealPoint Base = HLpoints[i, 0].Lerp(HLpoints[i, deviations + 1], (1 / (deviations+1)) * p);
                    double GradTheta = Math.Atan(HLpoints[i, 0].Grad(HLpoints[i, deviations + 1])); // +/- 0.5PI 

                    double swerve = R.NextDouble() * MaxDeviationSize;

                    HLpoints[i, p + 1] = Base;//new RealPoint(Base.x + (swerve * Math.Cos(GradTheta + ((0.5 * Math.PI) * ((R.Next(2) == 0) ? -1 : 1)))), Base.y + (swerve * Math.Sin(GradTheta + ((0.5 * Math.PI) * ((R.Next(2) == 0) ? -1 : 1)))));
                }
            }
        }

        public void Draw(Cursor c)
        {
            for (int i = 0; i < Hyperlanes; i++)
            {
                c.Line((float)HLpoints[i, 0].x, (float)HLpoints[i, 0].y, (float)HLpoints[i, deviations + 1].x, (float)HLpoints[i, deviations + 1].y, Color.Red);
                for (int p = 0; p < deviations+2; p++)
                {
                    c.Dot((float)HLpoints[i, p].x, (float)HLpoints[i, p].y, 2, 2, Color.Blue);
                }
            }
        }
    }

	public class Sector : FilledArc
	{
		public SubSector[,] SubSectors;

		public int SystemMult;
		public int[] TypeMult;

		public Sector(int tradius, int start, int end) : base(tradius, start, end) { }

		public void initSubsectors()
		{
			Debug.check();
			SubSectors = new SubSector[Consts.IsubDivisions,(End-Start)];
			for (int i = 0; i < Consts.IsubDivisions; i++)
			{
				for (int o = 0; o < (End - Start); o++)
				{
					SubSectors[i, o] = new SubSector(this, i, o+Start, (o + 1)+Start);
				}
            }
			Debug.p("Subsectors initialised in {0}\n", Debug.end());
		}

		public void Draw(Cursor c, bool DrawSubSectors = true)
		{
			Debug.check();
			if (DrawSubSectors)
			{
				for (int i = 0; i < SubSectors.GetLength(0); i++)
				{
					for (int o = 0; o < SubSectors.GetLength(1); o++)
					{
						SubSectors[i, o].Draw(c);
					}
				}
				Debug.p("Subsectors drawn in {0}\n", Debug.end());
				Debug.check();
			}

			c.Line(
				Consts.Radius + (Tradius * Consts.Division) * Math.Cos(Start * Consts.Segment),
				Consts.Radius + (Tradius * Consts.Division) * Math.Sin(Start * Consts.Segment),
				Consts.Radius + ((Tradius - 1) * Consts.Division) * Math.Cos(Start * Consts.Segment),
				Consts.Radius + ((Tradius - 1) * Consts.Division) * Math.Sin(Start * Consts.Segment)
            );
			c.Line(
				Consts.Radius + (Tradius * Consts.Division) * Math.Cos(End * Consts.Segment),
				Consts.Radius + (Tradius * Consts.Division) * Math.Sin(End * Consts.Segment),
				Consts.Radius + ((Tradius - 1) * Consts.Division) * Math.Cos(End * Consts.Segment),
				Consts.Radius + ((Tradius - 1) * Consts.Division) * Math.Sin(End * Consts.Segment)
			); //draw the sidelines

			c.Arc(Tradius, Start, End);
			c.Arc(Tradius-1, Start, End);
			Debug.p("Main outline drawn in {0}\n", Debug.end());
		}
	}

	public class SubSector : FilledArc
	{
		public Sector Parent;
        
		public int Type = 1; //0==deserted, 1==Normal, 2==popular, 3==nebula

		public SubSector(Sector p) { Parent = p; }
		public SubSector(Sector p, int tradius, int start, int end) : base(tradius,start,end) { Parent = p; }

		public void Draw(Cursor c)
		{
			double ThisTop = (Parent.Tradius * Consts.Division) - (Tradius * Consts.IDivision);
			double ThisBottom = (Parent.Tradius * Consts.Division) - ((Tradius + 1) * Consts.IDivision);

			c.Line(
				Consts.Radius + ThisTop * Math.Cos(Start * Consts.Segment),
				Consts.Radius + ThisTop * Math.Sin(Start * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Cos(Start * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Sin(Start * Consts.Segment),
				Color.CornflowerBlue
			);
			c.Line(
				Consts.Radius + ThisTop * Math.Cos(End * Consts.Segment),
				Consts.Radius + ThisTop * Math.Sin(End * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Cos(End * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Sin(End * Consts.Segment),
                Color.CornflowerBlue
			); //draw the sidelines

			double precision = 0.05;
			for (double i = Start; i < End; i += precision)
			{
				c.Line(
				Consts.Radius + ThisTop * Math.Cos(i * Consts.Segment),
				Consts.Radius + ThisTop * Math.Sin(i * Consts.Segment),
				Consts.Radius + ThisTop * Math.Cos((i + precision) * Consts.Segment),
				Consts.Radius + ThisTop * Math.Sin((i + precision) * Consts.Segment),
				Color.CornflowerBlue
				);
			}

			for (double i = Start; i <= End; i += precision)
			{
				c.Line(
				Consts.Radius + ThisBottom * Math.Cos(i * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Sin(i * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Cos((i + precision) * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Sin((i + precision) * Consts.Segment),
				Color.CornflowerBlue
				);
			}

			switch (Type)
			{
				case 0: Fill(c,Color.CornflowerBlue); break;

				case 2: Fill(c, Color.GreenYellow); break;

				case 3: Fill(c, Color.Purple); break;
			}
		}

		public void Fill(Cursor c, Color C)
		{
			double ThisTop = (Parent.Tradius * Consts.Division) - (Tradius * Consts.IDivision);
			double ThisBottom = (Parent.Tradius * Consts.Division) - ((Tradius + 1) * Consts.IDivision);

			double precision = 0.05;
			for (double i = Start; i < End; i += precision)
			{
				c.Line(
				Consts.Radius + ThisTop * Math.Cos(i * Consts.Segment),
				Consts.Radius + ThisTop * Math.Sin(i * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Cos(i * Consts.Segment),
				Consts.Radius + ThisBottom * Math.Sin(i * Consts.Segment),
				C
				);
			}
		}
	}

	public class Builder
	{
		Random r;

		int MinLength = 16; //minimum length of a sector in major divisions
		int MaxLength = 32; //maximum length of a sector in major divisions
		int CentreSize = 5; //The size of the empty space in the centre of the galaxy, in divisions

		int MinBiteLength = 2; //minimum bite length in Isubs
		int MaxBiteLength = 16; //maximum bite length in Isubs
		int MaxBiteHeight = 8; //maximum bite height in Hsubs

        Galaxy g;

		public Builder(int seed)
		{
			r = new Random(seed);
		}

		public Builder(string seed, Cursor c)
		{
			string TrueSeed = "";
			foreach (char i in seed)
			{
				TrueSeed += ((int)i*69).ToString();
			}
			Console.WriteLine(TrueSeed);
            r = new Random(TrueSeed.GetHashCode());
            g = new Galaxy(r);
            g.BuildHyperlanes();
            g.Draw(c);
		}

		public Builder()
		{
			r = new Random();
		}


		public Sector[] BuildMainSectors()
		{
			List<Sector> TempList = new List<Sector>();

			Debug.check();

			for (int i = Consts.Divisions; i > CentreSize; i--)
			{
				int Shift = r.Next(Consts.HsubDivisions);
				for (int o = Shift; o < Shift+Consts.HsubDivisions;)
				{
					int Length = ( o+MaxLength > Shift + Consts.HsubDivisions) ? (Shift+Consts.HsubDivisions)-o : r.Next(MinLength, MaxLength);
                    TempList.Add(new Sector(i,o,o+=Length));
				}
				Debug.p("Completed division {0} of {1} in {2}\n", Consts.Divisions - i, Consts.Divisions - CentreSize, Debug.reset());
			}
			Debug.end();
			return TempList.ToArray();
		}
		
		public void Bites(Sector Target, int Min, int Max)
		{
			for (int i = 0; i < r.Next(Min,Max)+1; i++)
			{
				Bite(Target);
			}
		}

		private void Bite(Sector ToBreak)
		{
			int Height = r.Next(1, MaxBiteHeight+1);
			int Length = r.Next(MinBiteLength, MaxBiteLength);

			int ActualLength = ToBreak.End - ToBreak.Start;

			int StartX = r.Next(0, ActualLength - Length);

			bool TB = r.Next(0, 2) == 1;

			for (int y = 0; y < Height; y++)
			{
				for (int x = 0; x < Length; x++)
				{
					ToBreak.SubSectors[y + ((TB) ? 0 : Consts.IsubDivisions - Height), StartX + x].Type = 0;
				}
			}
		}

		public void AssignRoles(ref Sector[] ToChange)
		{
		}
	}
*/
}
