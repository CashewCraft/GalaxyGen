using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace GalaxyGen
{
	public class Consts
	{
		//Begining of user defined inputs
		public const int Radius = 525;

		public const int Divisions = 16;

		public const int HsubDivisions = 512;
		public const int IsubDivisions = 16;
		
		//End of user defined inputs

		public const double Tau = Math.PI * 2;

		public const double Division = (Radius / Divisions);

		public const double Segment = (Tau / HsubDivisions);

		public const double IDivision = (Division/IsubDivisions);
	}

	public class Debug
	{
		public const bool Active = true;

		private static Stack<int> timestamps = new Stack<int>();

		public static void check()
		{
			if (Active)
			{
				timestamps.Push(Environment.TickCount);
			}
		}

		public static string reset()
		{
			if (Active)
			{
				int temp = timestamps.Pop();
				timestamps.Push(Environment.TickCount);
				return (Environment.TickCount - temp) + " ticks";
			}
		}

		public static string end()
		{
			if (Active)
			{
				return (Environment.TickCount - timestamps.Pop()) + " ticks";
			}
		}

		public static void p(string Word, params object[] args) { Console.Write(Word, args); }
	}

	public partial class Form1 : Form
	{
		public Form1()
		{
			InitializeComponent();
			this.Width = Consts.Radius*2;
			this.Height = Consts.Radius*2;
			this.BackColor = Color.Black;
		}

		private void Form1_Paint(object sender, PaintEventArgs e)
		{
			Cursor c = new Cursor(e, Color.Black);

			Builder b = new Builder("WH40K");

			Sector[] lis = b.BuildMainSectors();

			int count = 0;
			foreach (Sector a in lis)
			{
				Debug.check();
				a.initSubsectors();
				b.Bites(a,2,5);
				a.Draw(c);
				Debug.p("Drew sector #{0} in {1}\n\n",count++,Debug.end());
			}
		}
	}

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
		int Theta;
		int r;

		public Point(int R, int T)
		{
			Theta = T;
			r = R;
		}
	}

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

		Planet[] Planets;
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

	public class Planet
	{
		string Name;
		string StarType;
	}

	public class Builder
	{
		Random r;

		int MinLength = 16;
		int MaxLength = 32;
		int CentreSize = 5;

		int MinBiteLength = 2;
		int MaxBiteLength = 16;
		int MaxBiteHeight = 8;

		public Builder(int seed)
		{
			r = new Random(seed);
		}

		public Builder(string seed)
		{
			string TrueSeed = "";
			foreach (char i in seed)
			{
				TrueSeed += ((int)i*69).ToString();
			}
			Console.WriteLine(TrueSeed);
            r = new Random(TrueSeed.GetHashCode());
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
			foreach (Sector i in ToChange)
			{
				i.SystemMult = r.Next(1, 4);
				i.TypeMult = new int[0];//TODO: figure out who this is going to work
			}
		}
	}
}
