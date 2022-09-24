import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class findRetrovirus_RNAse_GAG_INT {
	static HashMap<String, String> codonTable = new HashMap<String, String>();
	static HashMap<String, String> origCandidates = new HashMap<String, String>();
	static HashMap<String, String[]> candidatesAmino = new HashMap<String, String[]>();
	static HashMap<String, String> candidates = new HashMap<String, String>();
	final static int NUM_RETRO_PATTERNS = 3;
	final static int NUM_RNA_PATTERNS = 2;
	final static int NUM_GAG_PATTERNS = 1;
	final static int NUM_INT_PATTERNS = 1;

	public static void main(String[] args) throws IOException {
		initCodonToAminoMap();
		loadChrMapTEMP(args);
		Pattern[] retroPatterns = new Pattern[NUM_RETRO_PATTERNS];
		Matcher[] retroMatchers = new Matcher[NUM_RETRO_PATTERNS];
		Pattern[] RNAPatterns = new Pattern[NUM_RNA_PATTERNS];
		Matcher[] RNAMatchers = new Matcher[NUM_RNA_PATTERNS];
		Pattern[] GAGPatterns = new Pattern[NUM_GAG_PATTERNS];
		Matcher[] GAGMatchers = new Matcher[NUM_GAG_PATTERNS];
		Pattern[] INTPatterns = new Pattern[NUM_INT_PATTERNS];
		Matcher[] INTMatchers = new Matcher[NUM_INT_PATTERNS];
		retroPatterns[0] = Pattern.compile("NAPA");
		retroPatterns[1] = Pattern.compile("KTAF");
		retroPatterns[2] = Pattern.compile("M.FGL");
		RNAPatterns[0] = Pattern.compile("DAS.........Q");
		RNAPatterns[1] = Pattern.compile("D..SR");
		GAGPatterns[0] = Pattern.compile(".............W..........................................L");
		INTPatterns[0] = Pattern.compile("H.......T");
		for (String key : candidatesAmino.keySet()) {
			for (String aminoSeq : candidatesAmino.get(key)) {
				retroMatchers[0] = retroPatterns[0].matcher(aminoSeq);
				retroMatchers[1] = retroPatterns[1].matcher(aminoSeq);
				retroMatchers[2] = retroPatterns[2].matcher(aminoSeq);
				RNAMatchers[0] = RNAPatterns[0].matcher(aminoSeq);
				RNAMatchers[1] = RNAPatterns[1].matcher(aminoSeq);
				GAGMatchers[0] = GAGPatterns[0].matcher(aminoSeq);
				INTMatchers[0] = INTPatterns[0].matcher(aminoSeq);
				for (Matcher retroMatcher : retroMatchers)
					if (retroMatcher.find()) {
						for(Matcher RNAMatcher : RNAMatchers) {
							if(RNAMatcher.find()) {
								for(Matcher GAGMatcher: GAGMatchers) {
									if(GAGMatcher.find()) {
										for(Matcher INTMatcher : INTMatchers) {
											if(INTMatcher.find()) {
												candidates.put(key, origCandidates.get(key));
											}
										}
									}
								}
							}
						}
					}
			}
		}
//		System.out.println(candidates.toString());
//		System.out.println(candidates.keySet().toString());
		writeNewCandidates(args);
	}

	private static void writeNewCandidates(String[] args) throws IOException {
		System.out.println(args[0]);
		String outputDirectory = args[0].replace("findPBS-sameleftright", "retrovirus-RNAase-GAG-INT-CandidateLTRs-");
		System.out.println(outputDirectory);
		File file = new File(outputDirectory);
		file.createNewFile();
		FileWriter writer = new FileWriter(file);
		for (String key : candidates.keySet()) 
			writer.write(key + "\n" + origCandidates.get(key) + "\n\n");
		
		writer.close();
	}

	public static void loadChrMapTEMP(String[] args) throws IOException {
		File chrs = new File(args[0]);
		Scanner fileReader = new Scanner(chrs);
		while (fileReader.hasNextLine()) {
			String name = fileReader.nextLine();
			String chr = fileReader.nextLine();
			fileReader.nextLine();
//			System.out.println(name + " -> " + chr);
			origCandidates.put(name, chr);
		}
		for (String key : origCandidates.keySet()) {
			candidatesAmino.put(key, convertBoth(origCandidates.get(key)));
//			System.out.println(key + " ---> \n" + Arrays.toString(convertBoth(candidates.get(key))));
		}
		fileReader.close();
	}

	public static void initCodonToAminoMap() {
		codonTable.put("TTT", "F");
		codonTable.put("ATT", "I");
		codonTable.put("GTT", "V");
		codonTable.put("TTC", "F");
		codonTable.put("CTC", "L");
		codonTable.put("ATC", "I");
		codonTable.put("GTC", "V");
		codonTable.put("TTA", "L");
		codonTable.put("CTA", "L");
		codonTable.put("ATA", "I");
		codonTable.put("GTA", "V");
		codonTable.put("TTG", "L");
		codonTable.put("CTG", "L");
		codonTable.put("ATG", "M");
		codonTable.put("GTG", "V");
		codonTable.put("TCT", "S");
		codonTable.put("CCT", "P");
		codonTable.put("ACT", "T");
		codonTable.put("GCT", "A");
		codonTable.put("TCC", "S");
		codonTable.put("CCC", "P");
		codonTable.put("ACC", "T");
		codonTable.put("GCC", "A");
		codonTable.put("TCA", "S");
		codonTable.put("CCA", "P");
		codonTable.put("ACA", "T");
		codonTable.put("GCA", "A");
		codonTable.put("TCG", "S");
		codonTable.put("CCG", "P");
		codonTable.put("ACG", "T");
		codonTable.put("GCG", "A");
		codonTable.put("TAT", "Y");
		codonTable.put("CAT", "H");
		codonTable.put("AAT", "N");
		codonTable.put("GAT", "D");
		codonTable.put("TAC", "Y");
		codonTable.put("CAC", "H");
		codonTable.put("AAC", "N");
		codonTable.put("GAC", "D");
		codonTable.put("CAA", "Q");
		codonTable.put("AAA", "K");
		codonTable.put("GAA", "E");
		codonTable.put("CAG", "Q");
		codonTable.put("AAG", "K");
		codonTable.put("GAG", "E");
		codonTable.put("TGT", "C");
		codonTable.put("CGT", "R");
		codonTable.put("AGT", "S");
		codonTable.put("GGT", "G");
		codonTable.put("TGC", "C");
		codonTable.put("CGC", "R");
		codonTable.put("AGC", "S");
		codonTable.put("GGC", "G");
		codonTable.put("CGA", "R");
		codonTable.put("AGA", "R");
		codonTable.put("GGA", "G");
		codonTable.put("TGG", "W");
		codonTable.put("CGG", "R");
		codonTable.put("AGG", "R");
		codonTable.put("GGG", "G");
		codonTable.put("CTT", "L");
		codonTable.put("TAA", "-");
		codonTable.put("TAG", "-");
		codonTable.put("TGA", "-");
	}

	public static String codonToAmino(String codon, boolean print) {
		if (print) {
//			System.out.println(codon + "->" + codonTable.get(codon));
//			String reversecodon = new StringBuilder(codon).reverse().toString();
//			System.out.println(reversecodon + "-->" + codonTable.get(reversecodon));
		}
		return codonTable.get(codon);
	}

	public static String codonToAmino(String codon) {
		return codonTable.get(codon);
	}

	public static String[] convertBoth(String codonSeq) {
		String[] out = new String[6];
		String genestrand1 = codonSeq;
		String genestrand2 = codonSeq.replace('A', 'H').replace('T', 'A').replace('H', 'T').replace('C', 'H')
				.replace('G', 'C').replace('H', 'G');
		String[] first = codonSeqToAminoSeq(genestrand1);
		String[] second = codonSeqToAminoSeqAntiPair(genestrand2);
		return new String[] { first[0], first[1], first[2], second[0], second[1], second[2] };
	}

	public static String[] codonSeqToAminoSeq(String codonSeq) {
		String aminoSeqF1 = "";
		String aminoSeqF2 = "";
		String aminoSeqF3 = "";
		for (int i = 0; i + 2 < codonSeq.length(); i += 3) {
			aminoSeqF1 += codonToAmino(codonSeq.substring(i, i + 3));
			if (i + 3 < codonSeq.length())
				aminoSeqF2 += codonToAmino(codonSeq.substring(i + 1, i + 4));
			if (i + 4 < codonSeq.length())
				aminoSeqF3 += codonToAmino(codonSeq.substring(i + 2, i + 5));
		}
		return new String[] { aminoSeqF1, aminoSeqF2, aminoSeqF3 };
	}

	public static String[] codonSeqToAminoSeqAntiPair(String codonSeq) {
		String aminoSeqF4 = "";
		String aminoSeqF5 = "";
		String aminoSeqF6 = "";
		codonSeq = new StringBuilder(codonSeq).reverse().toString();
		for (int i = 0; i + 2 < codonSeq.length(); i += 3) {
			aminoSeqF4 += codonToAmino(codonSeq.substring(i, i + 3));
			if (i + 3 < codonSeq.length())
				aminoSeqF5 += codonToAmino(codonSeq.substring(i + 1, i + 4));
			if (i + 4 < codonSeq.length())
				aminoSeqF6 += codonToAmino(codonSeq.substring(i + 2, i + 5));
		}

		return new String[] {aminoSeqF4, aminoSeqF5, aminoSeqF6};
	}
}
