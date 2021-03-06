package edu.stanford.facs.swing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class StochasticGradientDescent {
	public final int n_components;
	public double [][]head_embedding;
	public double [][]tail_embedding;
	private final int []head;
	private final int []tail;
	private final int n_epochs;
	private final int n_vertices;
	private final double []epochs_per_sample;
	private final double a;
	private final double b;
	private final double initial_alpha;
	public boolean move_other;
	public boolean []move_point;
	
	private double alpha;
	private final double BG2S;
	private final double ABNEG2;
	private final double BNEG1;
	private Random r;
	private final int n_1_simplices;
	private final double []epochs_per_negative_sample;
	private final double []epoch_of_next_negative_sample;
	private final double []epoch_of_next_sample;
	private double nTh;
	private int n_epoch;

	public static final boolean DEBUGGING=false;
	static boolean IsDebugIdx(int idx) {
		int []debug=new int[] {1654, 1707};
		for (int i=0;i<debug.length;i++) {
			if (debug[i]==idx) {
				return true;
			}
		}
		return false;
	}

	boolean doDebug =  false;
	int debugs=0;
	int computes=0;
	final int COMPUTE1=2, COMPUTE2=5;
	
	static void debug() {
		//StochasticGradientDescent sgd = getInstance(System.getProperty("user.home") + File.separator + "sgdInput.txt");
		StochasticGradientDescent sgd = getInstance(StochasticGradientDescent.class.getResourceAsStream("sgdInput.txt"));
		sgd.debugHeadTail1stLast();
		while (!sgd.nextEpochs()) {
		}
		System.out.println(sgd.n_components);
	}
	
	void debugHeadTail1stLast() {
		System.out.printf("%d/%d epochs:\thead 1st/end [%f %f]/[%f %f]\n",
				n_epoch, n_epochs,
				head_embedding[0][0],
				head_embedding[0][1],  
				head_embedding[head_embedding.length-1][0],
				head_embedding[head_embedding.length-1][1]);
		System.out.printf("\t\ttail 1st/end [%f %f]/[%f %f]\n",
				tail_embedding[0][0],
				tail_embedding[0][1],
				tail_embedding[tail_embedding.length-1][0],
				tail_embedding[tail_embedding.length-1][1]);
	}
	
	public static void main(String[] args) {
		
		if (true) {
			debug();
			return;
		}
		
		Double[][] answer = null;
	
		double []a= {5, 25, 4};
		double[]b=a;
		a[1]=19;
		System.out.println("a[1]="+a[1]+", b[1]="+b[1]);
		if(args.length>0) {
			answer = readfile(args[0]);
		}
		else { 
			answer = readfile("D:\\temp\\sample10k.csv");
		}
		System.out.println("There are " + answer.length + " rows");
		for(int row = 0; row < answer.length; row++) {
			for(int col = 0; col < answer[row].length; col++) {
				System.out.print(answer[row][col]);
				System.out.print(" ");
				
			}
			System.out.println();
			if(row > 10) {
				break;
			}
		}
	}
	
	public static StochasticGradientDescent getInstance(String filePath) {
		
		File file = new File(filePath);
		if (!file.exists() || !file.isFile()) {
			System.out.println(filePath + " does not exist or not a file");
			return null;
		}
		try {
			return getInstance(new FileInputStream(filePath));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return null;
		
	}
	
	public static StochasticGradientDescent getInstance(InputStream fr) {
	
		BufferedReader br = null;
        try {
    		br = new BufferedReader(new InputStreamReader(fr));
    		
            //head_embedding
    		int cols=Integer.parseInt(br.readLine());
            int rows = Integer.parseInt(br.readLine());
            
            double [][]head_embedding = new double[rows][cols];
            int arraySize = (rows*cols);
            for (int i=0; i < rows; i++) {
            	for (int j=0; j<cols;j++) {
            		head_embedding[i][j]=Double.parseDouble(br.readLine());
            	}
            }
            //tail_embedding
            rows = Integer.parseInt(br.readLine());
            arraySize = (rows*cols);
            double [][]tail_embedding = new double[rows][cols];
            for (int i=0; i < rows; i++) {
            	for (int j=0; j<cols;j++) {
            		tail_embedding[i][j]=Double.parseDouble(br.readLine());
            	}
            }
            //head
            rows = Integer.parseInt(br.readLine());
            int []head = new int[rows];
            for (int i=0; i < rows; i++) {
            	head[i]=Integer.parseInt(br.readLine());
            }
            //tail
            rows = Integer.parseInt(br.readLine());
            int []tail = new int[rows];
            for (int i=0; i < rows; i++) {
            	tail[i]=Integer.parseInt(br.readLine());
            }
            //n_epochs
            rows = Integer.parseInt(br.readLine());
            int n_epochs=Integer.parseInt(br.readLine());
            //n_vertices
            rows = Integer.parseInt(br.readLine());
            int n_vertices=Integer.parseInt(br.readLine());
            //epochs_per_sample
            rows = Integer.parseInt(br.readLine());
            double []epochs_per_sample = new double[rows];
            for (int i=0; i < rows; i++) {
            	epochs_per_sample[i]=Double.parseDouble(br.readLine());
            }
            //a
            rows = Integer.parseInt(br.readLine());
            double a=Double.parseDouble(br.readLine());
            //b
            rows = Integer.parseInt(br.readLine());
            double b=Double.parseDouble(br.readLine());
            //gamma
            rows = Integer.parseInt(br.readLine());
            double gamma=Double.parseDouble(br.readLine());
            //initial_alpha
            rows = Integer.parseInt(br.readLine());
            double initial_alpha=Double.parseDouble(br.readLine());
            //negative_sample_rate
            rows = Integer.parseInt(br.readLine());
            int negative_sample_rate=Integer.parseInt(br.readLine());
            StochasticGradientDescent sgd = new StochasticGradientDescent(head_embedding, 
            		tail_embedding, head, tail, n_epochs, n_vertices, epochs_per_sample, 
            		a, b, gamma, initial_alpha, negative_sample_rate);
            return sgd;
            
        } catch (final IOException e) {
        	e.printStackTrace();
        } finally {
        	if (br != null) {
        		try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
        	}
        }
    
		return null;
	}
	
	public static Double[][] readfile(String fileName){
		InputStreamReader ir = null;
		BufferedReader br = null;
        try
        {
        	if(fileName == null) { fileName = "D:\\temp\\sample10k.csv";}
        	File file = new File(fileName);
        	ir = new InputStreamReader(new FileInputStream(file), "UTF-8");
        	br = new BufferedReader(ir);

            // we don't know the amount of data ahead of time so we use lists
            List<Double> row = new ArrayList<>();
            List<Double[]> table = new ArrayList<>();
            final Double[][] none = new Double[0][];
            final Double[] none2 = new Double[0];
            
            String data = null;
            while ((data = br.readLine()) != null)
            {
                
                String [] arr = data.split(",");

                try {
                for(int i = 0; i<arr.length; i++) {
                	row.add(Double.parseDouble(arr[i]));
                }
                table.add(row.toArray(none2));
                row.clear();
                }
                catch(RuntimeException re) {
                	re.printStackTrace();
                }
            }
            return table.toArray(none);

        
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        finally
        {
            if (ir != null)
            {
                try {
					ir.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
            }
        }
        return null;
	}
	
	public void randomize() {
		r=new Random();
	}

	public double setReports(final double reports) {
		nTh=((double)n_epochs / reports);
		return nTh;
	}
	
	public double getReports() {
		return nTh;
	}
	
	public  StochasticGradientDescent(
			final double [][]head_embedding, final double [][]tail_embedding, 
			final int []head, final int []tail, final int n_epochs, final int n_vertices, 
			final double []epochs_per_sample, final double a, final double b, 
			final double gamma, final double initial_alpha, 
			final int negative_sample_rate){
		this.head_embedding=head_embedding;
		this.tail_embedding=tail_embedding;
		this.head=head;
		this.tail=tail;
		this.n_epochs=n_epochs;
		this.n_vertices=n_vertices;
		this.epochs_per_sample=epochs_per_sample;
		this.a=a;
		this.b=b;
		this.initial_alpha=initial_alpha;
		if (head_embedding.length >0) {
			this.n_components=head_embedding[0].length;
		}else {
			this.n_components=0;
		}
		move_other = head_embedding.length == tail_embedding.length;
		alpha = initial_alpha;
		BG2S=2*gamma* b;
		ABNEG2=-2.0*a*b;
		BNEG1=b-1;
		r=new Random(503l);
		n_1_simplices=epochs_per_sample.length;
		epochs_per_negative_sample=new double[epochs_per_sample.length];
		for (int i=0;i<n_1_simplices;i++) {
			epochs_per_negative_sample[i]=epochs_per_sample[i]/negative_sample_rate;
		}
		epoch_of_next_negative_sample=Arrays.copyOf(epochs_per_negative_sample, n_1_simplices);
		epoch_of_next_sample = Arrays.copyOf(epochs_per_sample, n_1_simplices);
		nTh=((double)n_epochs / (double)EPOCH_REPORTS());
		n_epoch=1;
	}
	
	public int getEpochsDone() {
		return n_epoch;
	}
	
	public int getEpochsToDo() {
		return n_epochs;
	}
	
	public boolean nextEpochs() {
		return nextEpochs(null);
	}
	boolean testJavaPointerToArray=true;
	
	public boolean nextEpochs( List<Integer>randis) {
		// start 1% speedup of making unchanging instance variables local final variables 
		if (DEBUGGING) {
			int []debugRands=null;
			debugRands= new int[]{6069,2754,4674,3957,446,2714,915,3983,4447,5739,598,8764};
			if (randis==null && debugRands != null) {
				randis=new ArrayList<>();
				for (int i=0;i<debugRands.length;i++) {
					randis.add(debugRands[i]);
				}
			}
		}
		final boolean move_other=this.move_other;
		final boolean []move_point=this.move_point;
		final int n_components=this.n_components;
		final int n_1_simplices=this.n_1_simplices;
		final int n_epochs=this.n_epochs;
		final int n_vertices=this.n_vertices;
		final int []head=this.head;
		final int []tail=this.tail;
		final double []epochs_per_sample=this.epochs_per_sample;
		final double a=this.a;
		final double b=this.b;
		final double initial_alpha=this.initial_alpha;
		final double BG2S=this.BG2S;
		final double ABNEG2=this.ABNEG2;
		final double BNEG1=this.BNEG1;
		final double []epochs_per_negative_sample=this.epochs_per_negative_sample;
		final double []epoch_of_next_negative_sample=this.epoch_of_next_negative_sample;
		final double []epoch_of_next_sample=this.epoch_of_next_sample;
		// finish 1% speedup of making final local variables
		
		int iRandi=0;
		double []current=new double[n_components];//= {0, 0};
		double []other=new double[n_components];//= {0, 0};
		int n_neg_samples=0;
		double []grad=new double[n_components];//={0, 0};
		double []sub=new double[n_components];//={0, 0};
		double grad_coef=0;
		double dist_squared=0;
		double val=0;
		if (move_other && move_point==null) { 
			//Ensure that head_embedding and tail_embedding point to the SAME memory
			//because identical matrices in MatLab are marshalled into JAVA as separate copies
			head_embedding=tail_embedding;
		}
		double alpha4=alpha*4, alphaNeg4=alpha*-4;
		
		
		for (int n=this.n_epoch;n<=n_epochs;n++) {
			for (int i=0;i<n_1_simplices;i++) {
				if (epoch_of_next_sample[i]>n) {
					continue;
				}
				final int j=head[i]-1;
				int k=tail[i]-1;
				for (int m=0;m<n_components;m++) {
					current[m]=head_embedding[j][m];
					other[m]=tail_embedding[k][m];
					sub[m]=current[m]-other[m];
				}
				dist_squared=0;
				for (int m=0;m<n_components;m++) {
					dist_squared+=sub[m]*sub[m];
				}
				if (dist_squared>0) {
					if (DEBUGGING) {
						//if (computes>=COMPUTE1 && computes<COMPUTE2 ){
						//doDebug=j==0 && n>7	 && n<10 && debugs<15;
						//doDebug=(j==1629) && n>5 && n<10  && debugs<9;
						//doDebug=(j==544) && n>5 && n<10  && debugs<9;
						//doDebug=(j==3983) && n>5 && n<10  && debugs<9;
						//doDebug=(j==3908 || j==3983) && (k==3908 || k==3983) && debugs<9;
						//doDebug=(j==2754) && n>5 && debugs<9;
						//doDebug=(j==8984) && n>5 && debugs<9;
						//doDebug=(j==3816) && n>5 && debugs<9;
						//doDebug=(j==1707) && n>5 && debugs<9;
						doDebug=(IsDebugIdx(j) || IsDebugIdx(k)) && n>5 && debugs<9;
						if (doDebug) {
							debugs++;
							System.out.printf("n=%d, j=%d, k=%d, sub=[%f %f], current=[%f %f], other=[%f %f], dist_squared=%f \n ", n, j, k, sub[0], sub[1], current[0], current[1], other[0], other[1], dist_squared);
						}
					}

					grad_coef=(ABNEG2*java.lang.Math.pow(dist_squared, BNEG1))/(a*java.lang.Math.pow(dist_squared, b)+1);
					for (int m=0;m<n_components;m++) {
						val=grad_coef*sub[m];
						if (val>=4) {
							grad[m]=alpha4;
						} else if (val <= -4) {
							grad[m]=alphaNeg4;
						} else {
							grad[m]=val*alpha;
						}
						current[m]=current[m]+grad[m];
					}
					if (DEBUGGING) {
						//if (computes>=COMPUTE1 && computes<COMPUTE2 ){
						if (doDebug) {
							System.out.printf(" ... grad=[%f %f], current=[%f %f]\n ", grad[0], grad[1], current[0], current[1]);
						}
					}
					if (move_other) {
						if ((move_point == null)) { 
							for (int m=0;m<n_components;m++) {
								other[m]=other[m]-grad[m];
								//Remove the next 2 comments to prove that head_embedding 
								//	and tail_embedding are pointers to the same memory
								//  and thus when one is changed so is the other
								
								//final double prior=head_embedding[k][m];
								
								tail_embedding[k][m]=other[m];//head_embedding references same memory
								
								/*
								if (testJavaPointerToArray) {
									if (prior != head_embedding[k][m]) {
										if (head_embedding[k][m]==tail_embedding[k][m]) {
											System.out.println("JAVA POINTER TO ARRAY TEST confirms understanding");
										}
									}
								}
								testJavaPointerToArray=false;
								*/
							}
						} else if (move_point[k]) {
							for (int m=0;m<n_components;m++) {
								other[m]=other[m]-grad[m];
								head_embedding[k][m]=other[m];
								tail_embedding[k][m]=other[m];
							}
						}
					}
				}else {
					for (int m=0;m<n_components;m++) {
						grad[m]=0;
					}
				}
				epoch_of_next_sample[i]+=epochs_per_sample[i];
				n_neg_samples = (int)Math.floor((((double)n) - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i]);
				for (int p=0;p<n_neg_samples;p++) {
					if (randis==null) {
						k=r.nextInt(n_vertices);
						//System.out.print(""+k+",");
					} else {
						if (iRandi>=randis.size())
							iRandi=0;
						k=randis.get(iRandi++);
					}
					if (move_other) {
						if ((move_point == null || move_point[k])) {
							if (j==k) {
								continue;
							}
						}
					}
					
					dist_squared=0;
					for (int m=0;m<n_components;m++) {
						other[m]=tail_embedding[k][m];
						sub[m]=current[m]-other[m];
						dist_squared+=sub[m]*sub[m];
					}
					if (dist_squared>0) {
						grad_coef=((BG2S/(0.001+dist_squared)))/(a*java.lang.Math.pow(dist_squared, b)+1);
						for (int m=0;m<n_components;m++) {
							val=grad_coef*sub[m];
							if (val>=4) {
								grad[m]=alpha4;
							} else if (val <= -4) {
								grad[m]=alphaNeg4;
							} else {
								grad[m]=val*alpha;
							}
						}
					} else {
						for (int m=0;m<n_components;m++) {
							grad[m]=4;
						}
					}
					for (int m=0;m<n_components;m++) {
						current[m]=current[m]+(grad[m]);
					}
					if (DEBUGGING) {
						if (doDebug){
							System.out.printf("\tp=%d, k=%d: dist_squared=%f, grad=[%f %f], \n\t\tcurrent=[%f %f]", 
									p, k, 
									dist_squared,
									grad[0],
									grad[1], 
									current[0],
									current[1]);
							System.out.print("\n");
						}
					}
				}
				for (int m=0;m<n_components;m++) {
					head_embedding[j][m]=current[m];
					if (move_point != null && move_point[j]) {
						tail_embedding[j][m]=current[m];
					}
				}
				//epoch_of_next_negative_sample(i) = epoch_of_next_negative_sample(i)+(n_neg_samples * epochs_per_negative_sample(i));
				epoch_of_next_negative_sample[i]+=n_neg_samples*epochs_per_negative_sample[i];
				if (DEBUGGING) {
					computes++;
				}
			}
			alpha = initial_alpha * (1 - (double)((double)n/(double)n_epochs));
			alpha4 = alpha*4;
			alphaNeg4=alpha*-4;
			/*
			if (DEBUGGING && n<23){
				System.out.printf("alpha %f, 4=%f, neg4=%f\n",alpha, alpha4, alphaNeg4);
			}
			*/
			if (Math.floor(((double)n)%nTh)==0) {
				this.n_epoch=n+1;
				if (DEBUGGING) {
					debugHeadTail1stLast();
					return true;
				}
				if (this.n_epoch<this.n_epochs) {
					return false;
				} else {
					return true;
				}
			}
		}
		return true;
	}
	
	public double [][]getEmbedding(){
		return head_embedding;
	}
	
	public boolean isFinished() {
		return this.n_epoch>=this.n_epochs;
	}
	
	final static boolean DEBUG_RANDI=false;
	final static boolean DEBUG_STATIC_AND_INSTANCE=false;
	public static final int EPOCH_REPORTS() {
		return 20;
	}
	
	public static double[][]Copy(final double[][]in){
		final int rows=in.length;
		final double[][]out=new double[in.length][];
		for (int row=0;row<rows;row++) {
			final int cols=in[row].length;
			out[row]=new double[cols];
			for (int col=0;col<cols;col++) {
				out[row][col]=in[row][col];
			}
		}
		return out;
	}
	
}
