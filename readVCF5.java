//no encoding/decoding stuff. so all doubles in kmer stuff turned to ints
package iproPackage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;

public class readVCF5 {
	//no changes to these 4
	public static final String REFGENOMEFILEPATH = "src/Mtub_H37Rv_NC000962.3.fasta";
	public static final int MTBREFLENGTH = 4411532;
	public static final String OUTPUT_FASTA_SUFFIX = ".fasta.txt";
    public static final String OUTPUT_KMER_SUFFIX = ".kmer.csv";
	//change  this every strain
    //public static final String VCFFILEPATH = "src/VCF/site.02.subj.0001.lab.2014222001.iso.1.vcf";
    //KMAX is the end of how many times we should get the kmer count for each vcf
	public static final int KMAX = 6;
	public static int c = 0;//muts w options
	public static int c2 = 0;//perfect muts
	public static int c3 = 0;//total muts
	public static int c4 = 0;//ignored muts
	//output paths are now created by the the program, to be named after the vcf file.
	//this code assumes the vcf directory created by August's python code has been put into the src folder.
	
    public static void main(String[] args) {
    	long startTime = System.currentTimeMillis();
    	char[] referenceGenome = new char[MTBREFLENGTH];
    	referenceGenome = readReferenceGenome(REFGENOMEFILEPATH, referenceGenome);
    	
    	String directoryPath = "C:\\Users\\crazy\\Downloads\\IPRO-main\\VCF";
    	
    	
    	try {
            // Get the directory path as a Path object
            Path directory = FileSystems.getDefault().getPath(directoryPath);

            // List all files in the directory
            List<Path> files = Files.list(directory).toList();

            int i = 0;
            // Iterate through the files
            for (Path file : files) {
                // Check if the file has a ".vcf" extension
                if (i == 5000 && i < 12287 && file.toString().endsWith(".vcf")) {
                    // Process each VCF file
                	System.out.printf("VCF File found: %3d%s\n",i, file.toString());
                	makeOutputs(referenceGenome, file.toString());
                }i++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    	
        
        //optional Console histogram data
        //histogramData(sortedKmers);
    	double percent = (double)c2 / (double)c;
        System.out.printf("Total # of Mutations: %d\nCommas: %3d\ngood splits: %3d\nPercent: %3f\n how many were missed? %d",c3,c,c2,percent,c4);
        //Console execution time
        long milliseconds = System.currentTimeMillis() - startTime;
        long minutes = milliseconds/60000;
        long seconds = milliseconds%60000 / 1000;
        System.out.println("Execution time: " + minutes + " minutes	" + seconds + " seconds");
    }

	private static void makeOutputs(char[] referenceGenome, String VCFFILEPATH) {
		int[] position;
        String[] ref, mut, mutA, mutB;
        
        //load vcf, get output names from vcf, print out output file names
        String vcfFileName = VCFFILEPATH.substring(39); // Remove "src/vcf/" prefix
        String outputFilePrefix = extractFileNamePrefix(vcfFileName);
        String directoryLocation = makeDirectory(outputFilePrefix);
        //String fastaFileA = directoryLocation+"/" + outputFilePrefix + ".fastaA.txt";        
        //String fastaFileB = directoryLocation+"/" + outputFilePrefix + ".fastaB.txt";

        //get the vcf file into an array, then split the array
        String[] vcfArray = readVCFFile(VCFFILEPATH);
        String[][] splitVCFArray = splitVCFStrings(vcfArray);

        //make the data types of size vcfArray.size
        position = new int[vcfArray.length];
        ref = new String[vcfArray.length];
        mut = new String[vcfArray.length];
        mutA = new String[vcfArray.length];
        mutB = new String[vcfArray.length];

        //process the data from splitVCFArray
        processVCFData(splitVCFArray, position, ref, mut);
        splitMutArray(mut, mutA, mutB);

        //find the total offset, then make the char arrays for the sequence data
        int totalOffsetA = calculateTotalOffset(ref, mutA);
        int totalOffsetB = calculateTotalOffset(ref, mutB);
        System.out.printf("Variation Offset Amount: A = %5d\t B = %5d\n", totalOffsetA, totalOffsetB);
        
        char[] strainA = new char[MTBREFLENGTH + totalOffsetA];
        char[] strainB = new char[MTBREFLENGTH + totalOffsetB];
        
        //load the data into the char arrays
        
        strainA = updatedReferenceGenome(strainA, referenceGenome, position, ref, mutA);
        strainB = updatedReferenceGenome(strainB, referenceGenome, position, ref, mutB);
        
        //choose which format for the output fasta file
        //writeCharArrayToFile(strainA, fastaFileA);
        //writeCharArrayToFile(strainA, fastaFileB);
        //writeCharArrayToGENBANKFile(strainA, fastaFile);
        
        //do kmer stuff
        for (int k = 1; k <= KMAX; k++) {
        	Map<String, Integer> kmerCounts = countKmers(strainA, k);
        	//List<Map.Entry<String, Integer>> sortedKmers = sortKmersByCount(kmerCounts);
        	List<Map.Entry<String, Integer>> sortedKmers = sortKmersAlphabetically(kmerCounts);
        	String kmerFileA = directoryLocation+"/" + outputFilePrefix + ".kmer"+k+"A.csv";
        	printKmerCountFile(sortedKmers, kmerFileA);
        	//do again for B
        	kmerCounts = countKmers(strainB, k);
        	//sortedKmers = sortKmersByCount(kmerCounts);
        	sortedKmers = sortKmersAlphabetically(kmerCounts);
        	String kmerFileB = directoryLocation+"/" + outputFilePrefix + ".kmer"+k+"B.csv";
        	printKmerCountFile(sortedKmers, kmerFileB);
        	//System.out.printf("Output1: %40s\nOutput2: %40s\n", fastaFileA, kmerFileA);
        }
	}

	public static String extractFileNamePrefix(String fileName) {
    	//takes in the name of vcf file, returns the prefix (uniqueID)
    	String ret = "";
    	String[] parts = fileName.split("\\."); //splits by '.' char
    	/*for (int i = 0; i < parts.length-2; i++) {
    		ret = String.join(".", ret, parts[i]);
    	}*/
    	for (int i = 0; i < parts.length-1; i++) {
    		ret = String.join(".", ret, parts[i]);
    	}
    	ret = ret.substring(1);
    	return ret;
    }
    
	private static String makeDirectory(String outputFilePrefix) {
		//makes the directories (if they do not exist already) names after uniqueID
    	String directoryLocation = "src/outputs/" + outputFilePrefix;
        File directory = new File(directoryLocation);
        if (!directory.exists()) {
            boolean created = directory.mkdirs();
            if (created) {
                System.out.println("Directory created successfully.");
            } else {
                System.out.println("Failed to create directory.");
            }
        } else {
            System.out.println("Directory already exists.");
        }
        return directoryLocation;
	}

	private static String[] readVCFFile(String filePath) {
		//takes in the filepath, returns an array with each line of the vcf file
        try {
            File file = new File(filePath);
            Scanner scanner = new Scanner(file);
            int vCount = 0;

            //first counts how many lines are in the vcf file
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                if (!line.contains("#")) { // Count non-commented lines
                    vCount++;
                }
            }
            scanner.close();

            //make the array the size of the vcf file
            String[] vcfArray = new String[vCount];
            scanner = new Scanner(file); // Reset scanner to start reading from the beginning

            //put vcf file data line by line into this array
            int nCount = 0;
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                if (!line.contains("#")) { // Skip commented lines
                    vcfArray[nCount++] = line;
                }
            }

            scanner.close();
            return vcfArray;
        } catch (FileNotFoundException e) {
            System.out.println("File not found: " + e.getMessage());
            return new String[0]; // Return empty array if file not found
        }
    }

    private static String[][] splitVCFStrings(String[] vcfArray) {
    	//takes in the array of each line from the vcf file, splits it into a 2d array
        String[][] splitVCFArray = new String[vcfArray.length][];
        for (int i = 0; i < vcfArray.length; i++) {
            splitVCFArray[i] = vcfArray[i].split("\t");
        }
        return splitVCFArray;
    }

    private static void processVCFData(String[][] splitVCFArray, int[] position, String[] ref, String[] mut) {
    	//processes the vcf file data from the split array. Takes in the POS REF and ALT columns
        for (int i = 0; i < splitVCFArray.length; i++) {
            position[i] = Integer.parseInt(splitVCFArray[i][1]);
            ref[i] = splitVCFArray[i][3];
            mut[i] = splitVCFArray[i][4];
            c3++;//total muts
            if (mut[i].contains(",")) {
            	//if the mutation has options, take the first 2 options
            	//System.out.printf("Comma here @line# %6d, mutation: %s\n", i+19, mut[i]);
            	c++; //num rows w commas
            	if (mut[i].split(",").length == 2)
            		c2++; //num perfect muts
            	else if (mut[i].split(",").length > 2)
            		c4 += mut[i].split(",").length - 2;//num missed muts 
            }
            //POS in vcf file is 1 indexed, change to zero indexed
            position[i] = position[i] - 1;
        }
        
    }
    
    private static void splitMutArray(String[] mut, String[] mutA, String[] mutB) {
		for (int i = 0; i < mut.length; i++) {
			if (mut[i].contains(",")) {
				mutA[i] = mut[i].split(",")[0];
            	mutB[i] = mut[i].split(",")[1];
			} else {
				mutA[i] = mut[i];
            	mutB[i] = mut[i];
			}
		}	
	}

	private static int calculateTotalOffset(String[] ref, String[] mut) {
		//compare each mutation to its reference for insertions and deletions
        int totalOffset = 0;
        for (int i = 0; i < ref.length; i++) {
        	//System.out.printf("line: %8d, mut.L:%5d, ref.L:%5d\n", i+19, mut[i].length(), ref[i].length());
        	totalOffset += mut[i].length();
        	totalOffset -= ref[i].length();
        }
        return totalOffset;
    }
    
    public static char[] readReferenceGenome(String filePath, char[] referenceGenome) {
    	//read in the reference genome fasta, store into a char array
        try {
        	int count = 0;
            File file = new File(filePath);
            Scanner scanner = new Scanner(file);

            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();

                if (line.equalsIgnoreCase(">NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome")) {
                	//check the first line is the right genome, then skip to the rest of the genome
                    System.out.println("Started Successfully");
                } else {
                	//copy each character over to the reference genome char array
                    for (int i = 0; i < line.length(); i++) {
                        referenceGenome[count + i] = line.charAt(i);
                    }
                    //count keeps track of each char that has already been counted to avoid overwriting
                    count += line.length();
                }
            }
            scanner.close();
        } catch (FileNotFoundException e) {
            System.out.println("File not found: " + e.getMessage());
        }
        return referenceGenome;
    }
    private static char[] updatedReferenceGenome(char[] strain, char[] referenceGenome, int[] position, String[] ref, String[] mut) {
    	//take in the reference genome, mutations and positions
    	//copy over reference to the strain genome if there is no mutation
    	//if there is a mutation, jump ahead the reference genome tracker by REF length
    	//and copy the mutation sequence to the strain genome
    	int i = 0; //strain
    	int j = 0; //referenceGenome
    	int k = 0; //vcf data
    	for (j = 0; j < referenceGenome.length; j++) {
    		if (j != position[k]) {//if no mutation
    			strain[i] = referenceGenome[j];
    			i++;
    		} else {
    			for (int x = 0; x < mut[k].length(); x++) {
    				strain[i] = mut[k].charAt(x);
    				i++;
    			}
    			j+=ref[k].length()-1;//skip over length of ref(should be same characters)
    			if (k<position.length-1)
    				k++;
    		}
    	}
    	//print out size of strain, size of reference, and number of mutations to console
    	System.out.printf("i = %8d, j = %8d, k = %8d\n", i,j,k);
        return strain;
    }
    
    private static void writeCharArrayToFile(char[] strain, String filePath) {
    	//print the strain sequence to the output file (STANDARD)
        try {
            FileWriter writer = new FileWriter(filePath);
            for(int i = 0; i < strain.length; i++) {
            	writer.write(strain[i]);
            	if (i%70 == 69)
            		writer.write('\n');
            	}
            
            writer.close();
            System.out.println("FASTA Seq Data written to file successfully");
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    }
    
    private static void writeCharArrayToGENBANKFile(char[] newFasta, String filePath) {
    	//print the strain sequence to the output file (GENBANK)
        try {
            PrintWriter writer = new PrintWriter(new FileWriter(filePath));
            for(int i = 0; i < newFasta.length; i++) {
            	writer.write(newFasta[i]);
            	if (i%10 == 9)
            		writer.write(' ');
            	if (i%60 == 59) {
            		writer.write('\n');
            		writer.printf("%10d ",i+2);
            	}
            }
            writer.close();
            System.out.println("FASTA Seq Data written to file successfully");
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    }
    public static Map<String, Integer> countKmers(char[] strain, int k) {//yikes this sucks
    	//count kmers based on constant int defined above main
    	//manually decode the special characters, if special character, only add .5 to the count
    	//special characters reference see: https://www.bioinformatics.org/sms/iupac.html
    	//this code sucks sorry
        Map<String, Integer> kmerCounts = new HashMap<>();

        for (int i = 0; i <= strain.length - k; i++) {
            String kmer = new String(strain, i, k);
            kmerCounts.put(kmer, kmerCounts.getOrDefault(kmer, 0) + 1);
        }
        return kmerCounts;
    }

	public static List<Map.Entry<String, Integer>> sortKmersByCount(Map<String, Integer> kmerCounts) {//chatgptmethod
        //sorts the kmers according to frequency
		List<Map.Entry<String, Integer>> sortedKmers = new ArrayList<>(kmerCounts.entrySet());
        sortedKmers.sort((a, b) -> b.getValue().compareTo(a.getValue()));
        return sortedKmers;
    }
	
	public static List<Map.Entry<String, Integer>> sortKmersAlphabetically(Map<String, Integer> kmerCounts) {
        List<Map.Entry<String, Integer>> sortedKmers = new ArrayList<>(kmerCounts.entrySet());

        // Comparator to sort alphabetically by key
        Comparator<Map.Entry<String, Integer>> comparator = Comparator.comparing(Map.Entry::getKey);

        sortedKmers.sort(comparator);
        return sortedKmers;
    }
	
	private static void printKmerCountFile(List<Entry<String, Integer>> sortedKmers, String filePath) {
		//print the sorted kmer count data to the output file
    	try {
            PrintWriter writer = new PrintWriter(new FileWriter(filePath));
            for (Map.Entry<String, Integer> entry : sortedKmers) {
                //System.out.println(entry.getKey() + ": " + entry.getValue());
            	if (entry.getValue() > 0) {
            		String row = entry.getKey() + "," + entry.getValue();
                    writer.println(row);
            	}
            }
            writer.close();
            //System.out.println("KMERCount Data written to file successfully");
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    	/*for (Map.Entry<String, Double> entry : sortedKmers) {
            System.out.println(entry.getKey() + ": " + entry.getValue());
        }*/
	}
	private static void histogramData(List<Entry<String, Double>> sortedKmers) {
    	//optional method for this kmer fount histogram
    	Map<Double, Integer> frequencyHistogram = new HashMap<>();
        double freq;
        for (Map.Entry<String, Double> entry : sortedKmers) {
            freq = entry.getValue();
            frequencyHistogram.put(freq, frequencyHistogram.getOrDefault(freq, 0) + 1);
        }
        // Print the histogram
        for (Map.Entry<Double, Integer> entry : frequencyHistogram.entrySet()) {
        	System.out.printf("Frequency: %f\t Count: %d\n", entry.getKey(), entry.getValue());
        }
	}
}