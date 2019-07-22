//--------------------------------------------------------------------------------
// File: Stats.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Tools
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Threading;
    using System.Threading.Tasks;

    /// <summary>
    /// Core statistical methods
    /// </summary>
    public static class Stats
    {
        public struct Point2
        {
            public double x { get; set; }
            public double y { get; set; }
        }
        
        // constants
        static double a1 = 0.254829592;
        static double a2 = -0.284496736;
        static double a3 = 1.421413741;
        static double a4 = -1.453152027;
        static double a5 = 1.061405429;
        static double p = 0.3275911;

        public static Random RNG = new Random();

        /// <summary>
        /// Create a list of integers from start to end (inclusive) increasing by incr
        /// </summary>
        /// <param name="start">Start value (inclusive)</param>
        /// <param name="end">End value (inclusive)</param>
        /// <param name="incr">Step value</param>
        public static List<int> Sequence(int start, int end, int incr)
        {
            if (incr <= 0)
            {
                throw new Exception(string.Format("Invalid sequence from {0} to {1} by {2} requested", start, end, incr));
            }

            var l = new List<int>();
            for (int i = start; i <= end; i += incr)
            {
                l.Add(i);
            }
            return l;
        }

        /// <summary>
        /// Create a list of integers from start to end (inclusive) increasing by incr and rounding as appropriate
        /// </summary>
        /// <param name="start">Start value (inclusive)</param>
        /// <param name="end">End value (inclusive)</param>
        /// <param name="incr">Step value</param>
        public static List<int> Sequence(double start, double end, double incr)
        {
            if (incr <= 0)
            {
                throw new Exception(string.Format("Invalid sequence from {0} to {1} by {2} requested", start, end, incr));
            }

            var l = new List<int>();
            for (double i = start; i <= end; i += incr)
            {
                l.Add((int)Math.Round(i));
            }
            return l;
        }

        /// <summary>
        /// Sample the specified count and max.
        /// </summary>
        /// <param name="count">Count.</param>
        /// <param name="max">Max.</param>
        public static int[] SampleWithReplacement(int count, int max)
        {
            var l = new int[count];
            for (int i = 0; i < count; i++)
            {
                l[i] = RNG.Next(max);
            }

            return l;
        }

        /// <summary>
        /// Samples without replacement.
        /// </summary>
        /// <returns>The without replacement.</returns>
        /// <param name="count">Count.</param>
        /// <param name="max">Max.</param>
        public static int[] SampleWithoutReplacement(int count, int max)
        {
            if (count > max)
            {
                throw new Exception(string.Format("Invalid sampling without replacement: count = {0} and max = {1}", count, max));
            }

            if (max > 10 * count)
            {
                HashSet<int> samples = new HashSet<int>();
                while (samples.Count < count)
                {
                    samples.Add(RNG.Next(max));
                }

                return samples.ToArray();
            }
            else
            {
                List<int> n = new List<int>();
                for (int i = 0; i < max; i++)
                {
                    n.Add(i);
                }

                if (count == max)
                {
                    int swap;
                    int index;
                    for (int i = 0; i < max; i++)
                    {
                        swap = n[i];
                        index = RNG.Next(max);
                        n[i] = n[index];
                        n[index] = swap;
                    }

                    return n.ToArray();
                }
                else
                {
                    int[] l = new int[count];
                    int index;
                    for (int i = 0; i < count; i++)
                    {
                        index = RNG.Next(n.Count);
                        l[i] = n[index];
                        n.RemoveAt(index);
                    }

                    return l;
                }
            }
        }

        /// <summary>
        /// Shuffle the specified count.
        /// </summary>
        /// <param name="count">Count.</param>
        public static int[] Shuffle(int count)
        {
            return SampleWithoutReplacement(count, count);
        }

        /// <summary>
        /// Sequences the array.
        /// </summary>
        /// <returns>The array.</returns>
        /// <param name="start">Start.</param>
        /// <param name="end">End.</param>
        public static int[] SequenceArray(int start, int end)
        {
            var a = new int[end - start + 1];

            for (int i = start; i <= end; i++)
            {
                a[i - start] = i;
            }

            return a;
        }

        public static double TwoTailedPvalue(double x)
        {
            double xa = Math.Abs(x);

            double pvalue = 2 * Phi(xa);

            if (xa > 6 || pvalue > 1 || pvalue < 0)
            {
                //Console.WriteLine(x + "\t" + xa + "\t" + pvalue);  
            }

            return pvalue;
        }

        public static double Phi(double x)
        {
            return 0.5 * errorfunctioncl(x);
        }

        public static double Log1Minus(double x)
        {
            double sum = 0;
            for (double n = 1; n < 1000.0; n++)
            {
                sum -= Math.Exp(Math.Log(x) * n - Math.Log(n));
            }

            return sum;
        }

        public static double PhiLowAccuracy(double x)
        {
            // Save the sign of x
            int sign = 1;
            if (x < 0)
                sign = -1;
            x = Math.Abs(x) / Math.Sqrt(2.0);
            
            // A&S formula 7.1.26
            double t = 1.0 / (1.0 + p*x);
            double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * Math.Exp(-x*x);

            return 0.5 * (1.0 + sign*y);
        }


        /// <summary>
        /// Log summation-based multiplication
        /// </summary>
        /// <returns>The m2.</returns>
        /// <param name="a">The alpha component.</param>
        /// <param name="b">The blue component.</param>
        /// <param name="x">x term of Phi function</param>
        private static double LSM2(double a, double b, double x)
        {
            return Math.Exp(Math.Log(a) + Math.Log(b) - (x * x));
        }

        /// <summary>
        /// Log summation-based multiplication of coefficient and base raised to power
        /// </summary>
        /// <returns>The p2.</returns>
        /// <param name="a">The alpha component.</param>
        /// <param name="b">The blue component.</param>
        /// <param name="pow">Pow.</param>
        private static double LSMP2(double a, double b, double pow, double x)
        {
            return Math.Exp(Math.Log(a) + pow * Math.Log(b) - (x * x));
        }

        /// <summary>
        /// Log subtraction-based division
        /// </summary>
        /// <returns>The m2.</returns>
        /// <param name="a">The alpha component.</param>
        /// <param name="b">The blue component.</param>
        private static double LBM2(double a, double b)
        {
            return Math.Exp(Math.Log(a) - Math.Log(b));
        }

        /// <summary>
        /// Log sum implementation of Phi function
        /// </summary>
        /// <returns>The L.</returns>
        /// <param name="x">The x coordinate.</param>
        public static double PhiLow(double x)
        {
            // Save the sign of x
            int sign = 1;
            if (x < 0)
                sign = -1;
            //x = Math.Abs(x) / Math.Sqrt(2.0);
            x = Math.Exp(Math.Log(Math.Abs(x)) - Math.Log(Math.Sqrt(2.0)));

            // A&S formula 7.1.26
            double t = Math.Exp(Math.Log(1.0) - Math.Log(1.0 + p*x));
            double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * Math.Exp(-x*x);


            //double y = 1.0 - ((((a5*t + a4)*t) + a3)*t + a2)*t*t + a1*t * Math.Exp(-x*x);
            //double y = 1.0 - (((a5*t + a4)*t) + a3)*t*t*t + a2*t*t + a1*t * Math.Exp(-x*x);
            //double y = 1.0 - ((a5*t + a4)*t)*t*t*t + a3*t*t*t + a2*t*t + a1*t * Math.Exp(-x*x);
            //double y = 1.0 - (a5*t*t + a4*t)*t*t*t + a3*t*t*t + a2*t*t + a1*t * Math.Exp(-x*x);
            double y2 = 1.0 - (a5*t*t*t*t*t * Math.Exp(-x*x) + a4*t*t*t*t * Math.Exp(-x*x) + a3*t*t*t * Math.Exp(-x*x) + a2*t*t * Math.Exp(-x*x) + a1*t * Math.Exp(-x*x));
            double y4 = 1.0 - (a5*t*t*t*t*t * Math.Exp(-x*x) + a4*t*t*t*t * Math.Exp(-x*x) + a3*t*t*t * Math.Exp(-x*x) + a2*t*t * Math.Exp(-x*x) + a1*t * Math.Exp(-x*x));

            double y3 = 1.0 - (LSMP2(a5, t, 5, x) - LSMP2(-a4, t, 4, x) + LSMP2(a3, t, 3, x) - LSMP2(-a2, t, 2, x) + LSM2(a1, t, x));

            if (0 == 1)
            {
                Console.WriteLine(y + "\t" + y2 + "\t" + y3);
                Console.WriteLine(y - y2);
                Console.WriteLine(y - y3);
                Console.WriteLine();
            };

            var pval =  0.5 * (1.0 + sign*y);

            if (1 == 0)
            {
                Console.WriteLine(y + "\t" + x + "\t" + y3);
            }

            return pval;
        }


        /*************************************************************************
        From alglib: Error function

        The integral is

                                  x
                                   -
                        2         | |          2
          erf(x)  =  --------     |    exp( - t  ) dt.
                     sqrt(pi)   | |
                                 -
                                  0

        For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
        erf(x) = 1 - erfc(x).


        ACCURACY:

                             Relative error:
        arithmetic   domain     # trials      peak         rms
           IEEE      0,1         30000       3.7e-16     1.0e-16

        Cephes Math Library Release 2.8:  June, 2000
        Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
        *************************************************************************/
        public static double errorfunction(double x)
        {
            double result = 0;
            double xsq = 0;
            double s = 0;
            double p = 0;
            double q = 0;

            s = Math.Sign(x);
            x = Math.Abs(x);
            if( (double)(x)<(double)(0.5) )
            {
                xsq = x*x;
                p = 0.007547728033418631287834;
                p = -0.288805137207594084924010+xsq*p;
                p = 14.3383842191748205576712+xsq*p;
                p = 38.0140318123903008244444+xsq*p;
                p = 3017.82788536507577809226+xsq*p;
                p = 7404.07142710151470082064+xsq*p;
                p = 80437.3630960840172832162+xsq*p;
                q = 0.0;
                q = 1.00000000000000000000000+xsq*q;
                q = 38.0190713951939403753468+xsq*q;
                q = 658.070155459240506326937+xsq*q;
                q = 6379.60017324428279487120+xsq*q;
                q = 34216.5257924628539769006+xsq*q;
                q = 80437.3630960840172826266+xsq*q;
                result = s*1.1283791670955125738961589031*x*p/q;
                return result;
            }
            if( (double)(x)>=(double)(10) )
            {
                result = s;
                return result;
            }
            result = s*(1-errorfunctionc(x));
            return result;
        }


        /*************************************************************************
        From alglib: Complementary error function

         1 - erf(x) =

                                  inf.
                                    -
                         2         | |          2
          erfc(x)  =  --------     |    exp( - t  ) dt
                      sqrt(pi)   | |
                                  -
                                   x


        For small x, erfc(x) = 1 - erf(x); otherwise rational
        approximations are computed.


        ACCURACY:

                             Relative error:
        arithmetic   domain     # trials      peak         rms
           IEEE      0,26.6417   30000       5.7e-14     1.5e-14

        Cephes Math Library Release 2.8:  June, 2000
        Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
        *************************************************************************/
        public static double errorfunctionc(double xd)
        {
            decimal x = (decimal)xd;

            double result = 0;
            decimal p = 0;
            decimal q = 0;

            if( (double)(x)<(double)(0) )
            {
                result = 2-errorfunctionc(-xd);
                return result;
            }
            if( (double)(x)<(double)(0.5) )
            {
                result = 1.0-errorfunction(xd);
                return result;
            }
            if( (double)(x)>=(double)(10) )
            {
                result = 0;
                return result;
            }
            p = 0.0m;
            p = 0.5641877825507397413087057563m+x*p;
            p = 9.675807882987265400604202961m+x*p;
            p = 77.08161730368428609781633646m+x*p;
            p = 368.5196154710010637133875746m+x*p;
            p = 1143.262070703886173606073338m+x*p;
            p = 2320.439590251635247384768711m+x*p;
            p = 2898.0293292167655611275846m+x*p;
            p = 1826.3348842295112592168999m+x*p;
            q = 1.0m;
            q = 17.14980943627607849376131193m+x*q;
            q = 137.1255960500622202878443578m+x*q;
            q = 661.7361207107653469211984771m+x*q;
            q = 2094.384367789539593790281779m+x*q;
            q = 4429.612803883682726711528526m+x*q;
            q = 6089.5424232724435504633068m+x*q;
            q = 4958.82756472114071495438422m+x*q;
            q = 1826.3348842295112595576438m+x*q;
            result = Math.Exp(-Math.Sqrt((double)x) + Math.Log((double)p) - Math.Log((double)q));
            return result;
        }

        private static double lsa(double a, double b)
        {
            return Math.Exp(Math.Log(a) + Math.Log(b));
        }

        public static double errorfunctioncl(double x)
        {
            double result = 0;
            double p = 0;
            double q = 0;

            if( (double)(x)<(double)(0) )
            {
                result = 2-errorfunctionc(-x);
                return result;
            }
            if( (double)(x)<(double)(0.5) )
            {
                result = 1.0-errorfunction(x);
                return result;
            }
            if( (double)(x)>=(double)(10) )
            {
                result = 0;
                return result;
            }
            p = 0.0;
            p = 0.5641877825507397413087057563+lsa(x,p);
            p = 9.675807882987265400604202961+lsa(x,p);
            p = 77.08161730368428609781633646+lsa(x,p);
            p = 368.5196154710010637133875746+lsa(x,p);
            p = 1143.262070703886173606073338+lsa(x,p);
            p = 2320.439590251635247384768711+lsa(x,p);
            p = 2898.0293292167655611275846+lsa(x,p);
            p = 1826.3348842295112592168999+lsa(x,p);
            q = 1.0;
            q = 17.14980943627607849376131193+lsa(x,q);
            q = 137.1255960500622202878443578+lsa(x,q);
            q = 661.7361207107653469211984771+lsa(x,q);
            q = 2094.384367789539593790281779+lsa(x,q);
            q = 4429.612803883682726711528526+lsa(x,q);
            q = 6089.5424232724435504633068+lsa(x,q);
            q = 4958.82756472114071495438422+lsa(x,q);
            q = 1826.3348842295112595576438+lsa(x,q);
            result = Math.Exp(-Math.Sqrt(x) + Math.Log(p) - Math.Log(q));
            return result;
        }

        public static double Pearson(IEnumerable<Point2> vector)
        {
            double Xbar = vector.Average(d => d.x);
            double Ybar = vector.Average(d => d.y);
            double productMoment = vector.Sum(d => (d.x - Xbar) * (d.y - Ybar));
            double devProduct = Math.Sqrt(vector.Sum(d => Math.Pow(d.x - Xbar, 2)) * 
                                          vector.Sum(d => Math.Pow(d.y - Ybar, 2)));
            
            return productMoment / devProduct;
        }

        public static double Pearson(List<double> x, List<double> y)
        {
            if (x.Count != y.Count)
            {
                throw new Exception("Invalid x and y vectors for pearson correlation; count mismatch");
            }

            double Xbar = x.Average();
            double Ybar = y.Average();
            double productMoment = x.Select((xi, i) => (xi - Xbar) * (y[i] - Ybar)).Sum();
            double devProduct = Math.Sqrt(x.Sum(xi => Math.Pow(xi - Xbar, 2)) * 
                y.Sum(yi => Math.Pow(yi - Ybar, 2)));

            return productMoment / devProduct;
        }
        
        public static double StudentTStatistic(double r, int n)
        {
            return r * Math.Sqrt( (double)(n - 2) / ( 1 - r * r) );
        }
        
        public static double FisherTransformZScore(double r, int n)
        {
            if (r == 1)
            {
                return double.MaxValue;
            }
            
            if (n <= 3)
            {
                throw new Exception("Invalid length for correlation");
            }
            
            return 0.5 * Math.Log( ( 1 + r ) / ( 1 - r ) ) * Math.Sqrt( n - 3 );
        }
        
        static double  NonparametricTest(int[] Tails, int[] SequenceLengths, Random rng)
        {
            int nSamples = 10000;
            
            int nObservedLessThanExpected = 0;
            
            for (int j = 0; j < nSamples; j++)
            {
                int[] randData = MultiSample(SequenceLengths, rng);
                int[] Differences = MultiDifference(Tails, randData).ToArray();
                int[] ResampledDifferences = MultiDifference(MultiSample(SequenceLengths, rng), randData).ToArray();
                
                if (Differences.Count(d => d < 0) - ResampledDifferences.Count(d => d < 0) < 0)
                {
                    nObservedLessThanExpected++;
                }
            }
            return 1.0 - (double)nObservedLessThanExpected / (double)nSamples;
        }
        
        public static IEnumerable<int> Shuffle(IEnumerable<int> data, Random rng)
        {
            int randRange = data.Count() * data.Count();
            return (from d in data
                    select new { Shuffle = rng.Next(randRange), Value = d }).OrderBy(x => x.Shuffle).Select(x => x.Value);
        }
        
        public static int[] MultiSample(IEnumerable<int> ranges, Random rng)
        {
            return (from r in ranges select rng.Next (r)).ToArray ();
        }

        /***************************************************************************************************************************************
         * Sampling and shuffling methods
         ***************************************************************************************************************************************/

        /// <summary>
        /// Sample with replacement
        /// </summary>
        /// <returns>The sample.</returns>
        /// <param name="n">N.</param>
        /// <param name="max">Max.</param>
        public static int[] BootstrapSample(int n, int max)
        {
            int[] s = new int[n];
            for (int i = 0; i < n; i++)
            {
                s[i] = RNG.Next(max);
            }
            
            return s;
        }

        /// <summary>
        /// Sample n elements from max (exclusive) without replacement
        /// </summary>
        /// <param name="n">number of samples</param>
        /// <param name="max">maximum value (exclusive)</param>
        public static int[] Sample(int n, int max)
        {
            var seq = Sequence(0, max - 1, 1);
            
            int[] samp = new int[n];
            for (int i = 0; i < n; i++)
            {
                // Sample from the first L elements
                int s = RNG.Next(max - i);
                
                // Extract sampled element and swap with tail element in sequence
                samp[i] = seq[s];
                seq[s] = seq[max - i - 1];
                seq[max - i - 1] = samp[i];
            }
            
            return samp;
        }

        /// <summary>
        /// Shuffle the integers from 1 to n (exclusive)
        /// </summary>
        /// <param name="n">Maximum integer (exclusive)</param>
        /*public static int[] Shuffle(int n)
        {
            int[] s = new int[n];
            for (int i = 0; i < n; i++)
            {
                s [i] = i;
            }
            
            for (int i = 0; i < n; i++)
            {
                int swap = s[i];
                int samp = RNG.Next(n);
                s[i] = s[samp]; 
                s[samp] = swap;
            }
            
            return s;
        }*.

        /***************************************************************************************************************************************
         * Bootstrapping methods
         ***************************************************************************************************************************************/

        /// <summary>
        /// Bootstrap the specified input data for n repetitions using the 
        /// stat function taking samples (with replacement) of size 
        /// </summary>
        /// <param name="input">Input.</param>
        /// <param name="n">N.</param>
        /// <param name="stat">Stat.</param>
        /// <param name="size">Size.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <typeparam name="Tret">The 2nd type parameter.</typeparam>
        public static Tret[] Bootstrap<T, Tret>(T[] input, int n, Func<T[], int[], Tret> stat, int size)
        {
            return BootstrapCore(input, n, stat, size, BootstrapSample);
        }

        /// <summary>
        /// Bootstrap the specified input data for n repetitions using the 
        /// stat function taking samples (without replacement) of size 
        /// </summary>
        /// <returns>The shuffle.</returns>
        /// <param name="input">Input.</param>
        /// <param name="n">N.</param>
        /// <param name="stat">Stat.</param>
        /// <param name="size">Size.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <typeparam name="Tret">The 2nd type parameter.</typeparam>
        public static Tret[] BootstrapShuffle<T, Tret>(T[] input, int n, Func<T[], int[], Tret> stat, int size)
        {
            return BootstrapCore(input, n, stat, size, (a, b) => Shuffle(a));
        }

        /// <summary>
        /// Parallels the bootstrap core allowing sampling to use different samplers
        /// </summary>
        /// <returns>The bootstrap core.</returns>
        /// <param name="input">Input.</param>
        /// <param name="n">N.</param>
        /// <param name="stat">Stat.</param>
        /// <param name="size">Size.</param>
        /// <param name="sampler">Sampler.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <typeparam name="Tret">The 2nd type parameter.</typeparam>
        private static Tret[] BootstrapCore<T, Tret>(T[] input, int n, Func<T[], int[], Tret> stat, int size, Func<int, int, int[]> sampler)
        {
            Tret[] data = new Tret[n];

            for (int i = 0; i < n; i++)
            {
                data[i] = stat(input, sampler(size, input.Length));
            }
            
            return data;
        }

        /***************************************************************************************************************************************
         * Parallel methods
         ***************************************************************************************************************************************/

        /// <summary>
        /// Parallel sampling, shuffling, and bootstrapping methods
        /// </summary>
        public class Threaded
        {
            /// <summary>
            /// Per-thread RNGs
            /// </summary>
            private static Dictionary<int, Random> threadRngs = new Dictionary<int, Random>();

            private static Queue<Random> unusedRngs = new Queue<Random>();

            /// <summary>
            /// Protected sampling rng accessor
            /// </summary>
            /// <returns>The rng.</returns>
            /// <param name="threadId">Thread identifier.</param>
            private static Random GetRng(int threadId)
            {
                lock (threadRngs)
                {
                    if (threadRngs.Count > 32)
                    {
                        unusedRngs = new Queue<Random>(threadRngs.Values);
                        threadRngs.Clear();
                    }

                    if (!threadRngs.ContainsKey(threadId))
                    {
                        Random rng = null;
                        if (unusedRngs.Count > 0)
                        {
                            rng = unusedRngs.Dequeue();
                        }
                        else
                        {
                            int seed = (int)(DateTime.Now.Ticks << 32);
                            seed = (seed << threadId) | (seed >> (32 - threadId));
                            rng = new Random(seed);
                        }

                        threadRngs.Add(threadId, rng);
                        //Console.WriteLine("{0} rngs used", threadRngs.Count);
                    }

                    return threadRngs[threadId];
                }
            }

            /***************************************************************************************************************************************
             * Parallel Sampling methods
             ***************************************************************************************************************************************/

            /// <summary>
            /// Sample with replacement
            /// </summary>
            /// <returns>The sample.</returns>
            /// <param name="n">N.</param>
            /// <param name="max">Max.</param>
            /// <param name="threadId"></param>
            public static int[] BootstrapSample(int n, int max, int threadId)
            {
                int[] s = new int[n];
                for (int i = 0; i < n; i++)
                {
                    s[i] = GetRng(Thread.CurrentThread.ManagedThreadId).Next(max);
                }
                
                return s;
            }
            
            /// <summary>
            /// Sample n elements from max (exclusive) without replacement
            /// </summary>
            /// <param name="n">number of samples</param>
            /// <param name="max">maximum value (exclusive)</param>
            /// 
            public static int[] Sample(int n, int max, int threadId)
            {
                var seq = Sequence(0, max - 1, 1);
                
                int[] samp = new int[n];
                for (int i = 0; i < n; i++)
                {
                    // Sample from the first L elements
                    int s = GetRng(Thread.CurrentThread.ManagedThreadId).Next(max - i);
                    
                    // Extract sampled element and swap with tail element in sequence
                    samp[i] = seq[s];
                    seq[s] = seq[max - i - 1];
                    seq[max - i - 1] = samp[i];
                }
                
                return samp;
            }
            
            /// <summary>
            /// Shuffle the integers from 1 to n (exclusive)
            /// </summary>
            /// <param name="n">Maximum integer (exclusive)</param>
            /// <param name="threadId"></param>
            public static int[] Shuffle(int n, int threadId)
            {
                //Console.WriteLine("Parallel shuffle in thread {0}, {1}", threadId, Thread.CurrentThread.ManagedThreadId);
                int[] s = new int[n];
                for (int i = 0; i < n; i++)
                {
                    s [i] = i;
                }
                
                for (int i = 0; i < n; i++)
                {
                    int swap = s[i];
                    int samp = GetRng(Thread.CurrentThread.ManagedThreadId).Next(n);
                    s[i] = s[samp]; 
                    s[samp] = swap;
                }
                
                return s;
            }

            /***************************************************************************************************************************************
             * Parallel Bootstrapping methods
             ***************************************************************************************************************************************/

            /// <summary>
            /// Parallel implementation to bootstrap the specified input data for n repetitions using the
            /// stat function taking samples (with replacement) of size 
            /// </summary>
            /// <param name="input">Input.</param>
            /// <param name="n">N.</param>
            /// <param name="stat">Stat.</param>
            /// <param name="size">Size.</param>
            /// <typeparam name="T">The 1st type parameter.</typeparam>
            /// <typeparam name="Tret">The 2nd type parameter.</typeparam>
            public static Tret[] Bootstrap<T, Tret>(T[] input, int n, Func<T[], int[], Tret> stat, int size)
            {
                return BootstrapCore(input, n, stat, size, BootstrapSample);
            }

            /// <summary>
            /// Parallel implementation to bootstrap the specified input data for n repetitions using the 
            /// stat function taking samples (without replacement) of size 
            /// </summary>
            /// <returns>The shuffle.</returns>
            /// <param name="input">Input.</param>
            /// <param name="n">N.</param>
            /// <param name="stat">Stat.</param>
            /// <param name="size">Size.</param>
            /// <typeparam name="T">The 1st type parameter.</typeparam>
            /// <typeparam name="Tret">The 2nd type parameter.</typeparam>
            public static Tret[] BootstrapShuffle<T, Tret>(T[] input, int n, Func<T[], int[], Tret> stat, int size)
            {
                return BootstrapCore(input, n, stat, size, (sampleSize, max, threadId) => Shuffle(sampleSize, threadId));
            }

            /// <summary>
            /// Parallels the bootstrap core allowing sampling to use different samplers
            /// </summary>
            /// <returns>The bootstrap core.</returns>
            /// <param name="input">Input.</param>
            /// <param name="n">N.</param>
            /// <param name="stat">Stat.</param>
            /// <param name="size">Size.</param>
            /// <param name="sampler">Sampler.</param>
            /// <typeparam name="T">The 1st type parameter.</typeparam>
            /// <typeparam name="Tret">The 2nd type parameter.</typeparam>
            private static Tret[] BootstrapCore<T, Tret>(T[] input, int n, Func<T[], int[], Tret> stat, int size, Func<int, int, int, int[]> sampler)
            {
                Tret[] data = new Tret[n];

                Parallel.For(0, n, i => 
                {
                    var id = Thread.CurrentThread.ManagedThreadId;
                    data[i] = stat(input, sampler(size, id, input.Length));
                    //stat(input, sampler(size, input.Length));
                });
                
                return data;
            }
        }


        
        public static IEnumerable<int> MultiDifference(int[] d1, int[] d2)
        {
            int[] d = new int[d1.Count()];
            if (d1.Count() != d2.Count())
            {
                throw new System.Exception("Entry sizes don't match!");
            }
            //Console.WriteLine("Size of d {0}", d1.Count);
            Parallel.For(0, d.Length, i => d[i] = d1[i] - d2[i]);
            //Console.WriteLine("\t{0}\t{1}\t{2}", d[0], d1[0], d2[0]);
            return d;
        }
        
        public class MotifInstance
        {
            public string Name { get; set; }
            public int Index { get; set; }
            public int Length { get; set; }
            public string Value { get; set; }
            public MotifInstance(){}
        }
        
        public static IEnumerable<Tuple<MotifInstance, MotifInstance>> OverlappingRanges(IEnumerable<MotifInstance> Matches1, IEnumerable<MotifInstance> Matches2)
        {
            return from m1 in Matches1
                from m2 in Matches2
                    where RangeOverlap(m1.Index, m1.Index + m1.Length - 1, m2.Index, m2.Index + m2.Length - 1)
                    select new Tuple<MotifInstance, MotifInstance>(m1, m2);
        }
        
        private static bool Between(int s, int e, int t)
        {
            return (t >= s && t <= e);
        }
        
        static bool RangeOverlap(int s1, int e1, int s2, int e2)
        {
            //Console.WriteLine("{0},{1},{2},{3}", s1, e1, s2, e2);
            return  (Between(s1, e1, s2) ||
                     Between(s1, e1, e2)) ||
                (Between(s2, e2, s1) ||
                 Between(s2, e2, e1));
        }

        /// <summary>
        /// Sample standard deviation
        /// </summary>
        /// <returns>The dev.</returns>
        /// <param name="data">Data.</param>
        public static double StdDev(IEnumerable<double> data)
        {
            if (data.Count() == 0)
            {
                return double.NaN;
            }

            if (data.Count() == 1)
            {
                return double.NaN;
            }

            double mean = data.Average();
            double ss = data.Sum(x => Math.Pow(x - mean, 2.0));
            return Math.Sqrt(ss / (data.Count() - 1));
        }

        /// <summary>
        /// Standard sample error
        /// </summary>
        /// <returns>The error.</returns>
        /// <param name="data">Data.</param>
        public static double StdErr(IEnumerable<double> data)
        {
            return StdDev(data) / Math.Sqrt(data.Count());
        }


        public static Tuple<double, double, double> Quartiles(IEnumerable<double> data)
        {
            var orderedData = data.OrderBy(x => x).ToList();
            int count = orderedData.Count;

            double median = Median(orderedData);
            double lower = Median(orderedData.Take(count / 2));
            double upper = Median(orderedData.Reverse<double>().Take(count / 2));

            return new Tuple<double, double, double>(lower, median, upper);
        }

        public static double Median(IEnumerable<double> data)
        {
            int count = data.Count();
            if (count == 0)
            {
                return double.NaN;
            }

            if (count % 2 == 1)
            {
                return data.OrderBy(x => x).Take(count / 2 + 1).Last();
            }
            else
            {
                return data.OrderBy(x => x).Take(count / 2 + 1).Reverse().Take(2).Average();
            }
        }

        public class StatisticalTests
        {
            public static double HypergeometricTest(int N, int n, int m, int k)
            {
                double flCdf = 0.0;
                if (k > n || k > m)
                {
                    throw new System.Exception(String.Format ("Hypergeometric overlap larger than either set size: n = {0}, m = {1}, k = {2}", n, m , k));
                }
                
                if (n > N || m > N)
                {
                    throw new System.Exception(String.Format ("Hypergeometric subsets larger than superset: n = {0}, m = {1}, N = {2}", n, m, N));
                }
                
                for (int i = k; i <= Math.Min(n, m); i++)
                {
                    flCdf += HypergeometricMass(N, n, m, i);
                }
                
                return flCdf;
            }
            
            public static double HypergeometricMass(int N, int n, int m, int k)
            {
                return Math.Exp(logHypergeometricMass(N, n, m, k));
            }
            
            public static double logHypergeometricMass(int N, int n, int m, int k)
            {
                return logBinomialCoef(m, k) + logBinomialCoef(N-m, n-k) - logBinomialCoef(N, n);
            }
            
            public static double BinomialCoef(int n, int k)
            {
                return Math.Exp(logBinomialCoef(n, k));
            }
            
            enum BinomialTerms { N, K, NmK };
            static double[] s_aflLogSums;
            
            public static double logBinomialCoef(int n, int k)
            {
                if (k > n)
                {
                    throw new System.Exception(String.Format ("Binomial coefficient with k > n: k = {0}, n = {1}", k, n));
                }
                
                if (s_aflLogSums == null || n >= s_aflLogSums.Length)
                {
                    double[] aflOldLogSums = s_aflLogSums;
                    s_aflLogSums = new double[2 * n];
                    if (aflOldLogSums != null)
                    {
                        Parallel.For(0, aflOldLogSums.Length, i => s_aflLogSums[i] = aflOldLogSums[i]);
                        double flLogSum = aflOldLogSums.Last();
                        for (int i = aflOldLogSums.Length; i < s_aflLogSums.Length; i++)
                        {
                            flLogSum += Math.Log (i);
                            s_aflLogSums[i] = flLogSum;
                        }
                    }
                    else
                    {
                        double flLogSum = 0.0;
                        for (int i = 1; i < s_aflLogSums.Length; i++)
                        {
                            flLogSum += Math.Log(i);
                            s_aflLogSums[i] = flLogSum;
                        }
                        
                    }
                }
                
                return s_aflLogSums[n] - (s_aflLogSums[k] + s_aflLogSums[n - k]);
            }


            public static double PermutationTest(double[] set1, double[] set2, int repCount)
            {
                double meanDiff = Math.Abs(set1.Average() - set2.Average());

                int l1 = set1.Length;
                int l2 = set2.Length;


                Func<double[], int[], int> test = (data, sample) => 
                {
                    double sum1 = 0;
                    for (int i = 0; i < l1; i++)
                    {
                        sum1 += data[sample[i]];
                    }
                    
                    double sum2 = 0;
                    for (int i = 0; i < l2; i++)
                    {
                        sum2 += data[sample[l2 + i]];
                    }
                    
                    return Math.Abs(sum1/l1 - sum2/l2) > meanDiff ? 1 : 0;
                };

                /*Func<double[], int[], int> t2 = (data, sample) =>
                {
                    return 1;
                };*/

                return (double)BootstrapShuffle(set1.Concat(set2).ToArray(), repCount, test, set1.Length + set2.Length).Sum() / (double)repCount;
            }
        }
    }
}

