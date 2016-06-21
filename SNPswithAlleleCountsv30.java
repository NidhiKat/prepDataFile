package prepDataFile;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;
/**
 * Update the attribute section as well as Data Section 
 * Error as omn 14th June : The data section is updated with only ?. 
 * @author skatla
 *
 */
public class SNPswithAlleleCountsv30 {
	// Name of the file to be written to 
		public final static String WRITE_TO_FILENAME = "snpAlleleCounts_June15v20_TreeMap.txt";
			
		/** Read Files from this location*/
		public final static String READ_FROM_PATH = "E:/SNPData/DataForAmblyopia/caseControl/";
		
		/** Final Arff file name*/
		public final static String PRINT_FINAL= "final_arff_15Junev20.arff";
		
		/** CALL RATE - SET at 95% of population*/
		public final static int CALL_RATE = (int)(0.95*68);
		
		/**Threshold value for MAF - SNPs with MAF less than threshold value will be removed */
		public final static double THRESHOLD_MAF = 0.01;
				
		/**  HashMap of RSID and alleles with their counts.*/
		public static HashMap<String, TreeMap<String, Integer>> snpAlleleCount = new HashMap<String, TreeMap<String, Integer>>();
		
		/** % missing value above which the subject's data is not used - This must be multipled by the snpAlleleCount size  */
		public static double THRESHOLD_MISSING = 0.05;
		
		/**  HashMap for rsid and SNP count*/
		public static HashMap<String, Double>SNPCount = new HashMap<String, Double>();
			
				
		/**
		 *  Outlier Detecting Method
		 *  Updated on 14 June	
		 * @param snpValue
		 * @param valuesMap
		 * @return
		 */
		public static boolean isOutlier(String snpValue, TreeMap<String, Integer> valuesMap){
			boolean result1 = true;
			boolean result2 = true;
			boolean result = true;
			if(snpValue.length() != 2){
				result = true;
				return result;
			}
			
			if(valuesMap.size() > 2){
				int maxVal = 0;	
				int minVal = 999999;
				for(String key :valuesMap.keySet()){
					if(!(key.equals(snpValue))){
						if(maxVal < valuesMap.get(key) )
							maxVal = valuesMap.get(key);
						if(minVal > valuesMap.get(key))
							minVal = valuesMap.get(key);
						if(result1 && (key.indexOf(snpValue.charAt(0)) != -1)){
							result1 = false;
						}
							
						if(result2 && (key.indexOf(snpValue.charAt(1)) != -1)){
							result2 = false;
						}
							
						if((!result1 )&& (!result2)){
							result = false;
							return result;
						}
							
					}
				}
				
				if((valuesMap.size() == 3) && (result1 || result2)){
					if((maxVal <= valuesMap.get(snpValue)) || (minVal < valuesMap.get(snpValue)))
						result = false;
				}
			}
			else 
				result = false;
			return result;			
		}
		
		
		/**
		 * Removes the outliers, the SNPs with zero or one Alleles, SNPs with call rate < 95 % (SNPs missing in more than 5% population)
		 * SNPCount holds the number of times the SNP had valid value. 
		 */
		public static void preProcessSNPAlleleCount(){
			Iterator<Map.Entry<String, TreeMap<String, Integer>>> it = snpAlleleCount.entrySet().iterator();
			while(it.hasNext()){
				Map.Entry<String, TreeMap<String, Integer>> ent = it.next();
				Iterator<Map.Entry<String, Integer>> iter = ent.getValue().entrySet().iterator();
				int valCount =0;
				while(iter.hasNext()){
					Map.Entry<String, Integer> entry = iter.next();
					if(isOutlier(entry.getKey(), ent.getValue()))
						iter.remove();
					else valCount += entry.getValue();
				}
				
				SNPCount.put(ent.getKey(), (double)valCount);
				
				if((ent.getValue().size() == 1) ||(ent.getValue().size() == 0) || (valCount < CALL_RATE) ){
					System.out.println("Remove Entry "+ent.getKey());
					it.remove();
				}
			}
		}
		
		
		
		/**
		 * For a given RSID, the alleles are adjusted to have format XX, XY,YY. 
		 * If one of them is missing, it is populated from the other two.
		 * Created on : 14 June
		 * @param valuesMap
		 * @param rsid
		 */
		
		public static void adjustSNPValues(TreeMap<String,Integer> valuesMap, String rsid){
			System.out.println("Enter adjustSNPValues" + rsid);
			String firstKey = valuesMap.firstKey();
			String secondKey = valuesMap.lastKey();
			System.out.println(firstKey);
			
		//	char firstChar = firstKey.charAt(0);
			char[] input = new char[2];
			if(firstKey.charAt(0) != firstKey.charAt(1)){
				input[0] = firstKey.charAt(0);
				input[1] = firstKey.charAt(0);
				
			} else if(secondKey.charAt(0) != secondKey.charAt(1)){
				input[0] = secondKey.charAt(1);
				input[1] = secondKey.charAt(1);
				
			} else{
				input[0] = firstKey.charAt(0);
				input[1] = secondKey.charAt(1);
			}
			valuesMap.put(String.valueOf(input), 0);
			snpAlleleCount.put(rsid, valuesMap);
		}
		
		
		
		/**
		 * Calculate the MAF value for each SNP. In main, we remove the SNPs that have MAF < 0.01
		 * Written on 14 June
		 * @param snpValueTree
		 * @param snpCount
		 * @return
		 */
		public static void calcMAF(TreeMap<String, Integer>snpValueTree, String rsid){
			double resultMAF = 0.0;
			double  minorAlleleCount = ((snpValueTree.get(snpValueTree.firstKey()) > snpValueTree.get(snpValueTree.lastKey()) ) ? snpValueTree.get(snpValueTree.lastKey()) : snpValueTree.get(snpValueTree.firstKey()) );
			minorAlleleCount = (minorAlleleCount*2)+(snpValueTree.higherEntry(snpValueTree.firstKey()).getValue());	
			try{
				 resultMAF = (double)(minorAlleleCount / (2*SNPCount.get(rsid)));
			//	 System.out.println(rsid+ "  "+minorAlleleCount +"  "+SNPCount.get(rsid));
			} catch(ArithmeticException e){
				System.out.println(e.getMessage());
				System.out.println("Error : SNPCount has 0 value at" + rsid );
			}
			SNPCount.put(rsid, resultMAF);
		//	return resultMAF;
		}
		
		/**
		 * Print the master Map along with SNP counts/ maf 
		 * @param args
		 */
		public static void printsnpAlleleCount(String filename){
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
				 snpAlleleCount.forEach((k,v) -> {
					 
					StringBuilder updLine = new StringBuilder("@attribute ");
					StringBuilder valueList = new StringBuilder();
					  	
					updLine.append(k).append(" {");
					
					v.forEach((k1,v1) -> {
						updLine.append("'").append(k1).append("',");
						valueList.append(v1).append(",");
					});
					
					updLine.append("} {").append(valueList).append("} ").append(SNPCount.get(k));
					try {
						bw.write(updLine.toString());
						bw.newLine();
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						e.getMessage();
						System.out.println("Error while writing to file");
					}
				 });
				 bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				e.getMessage();
				System.out.println(" Error while creating the file to write to ");
			}
		}
		
		
		
		
		public static void printSubjectsAndSNPs(String filename){
			try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			bw.write("@RELATION SNPs");
			bw.newLine();
			bw.write("@attribute user_ids { userIdLine");
			bw.newLine();
			
			 snpAlleleCount.forEach((k,v) -> {
				 
				StringBuilder updLine = new StringBuilder("@attribute ");
				StringBuilder valueList = new StringBuilder();
				  	
				updLine.append(k).append(" {");
				
				v.forEach((k1,v1) -> {
					updLine.append("'").append(k1).append("',");
					});
				
				updLine.append("}");
				try {
					bw.write(updLine.toString());
					bw.newLine();
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					e.getMessage();
					System.out.println("Error while writing to file");
				}
			 });
			 
			// Note that the "Attribute" is in capital letters here as against  small letter earlier and the code doesn't introduce an additional attribute while filling the user instance data. 
			// This requirement is this particular code specific and not weka specific. 
			bw.write("@ATTRIBUTE 'Class' {'Case', 'Control'}");
			bw.newLine();
			// Create the @DATA line
			bw.write("@DATA");
			bw.close();
			// Read each users' file in case(sets case-control to 1 and update it as an instance to allSNPsFinal.arff
			// send case before control 
			updateUserData("Case", "E:/SNPData/DataForAmblyopia/case");
			updateUserData("Control", "E:/SNPData/DataForAmblyopia/control");
			
		} catch(IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			e.getMessage();
			System.out.println(" Error while creating the file to write to ");
		}
		}
		
		public static void updateUserMap(HashMap<String, String> userMap){
			snpAlleleCount.forEach((k,v)->{
				userMap.put(k, "?");
			});
		}
		
		public static int countMissingValues(HashMap<String, String> rsidGenotype){
			int count_missing =0;
			for(String key : rsidGenotype.keySet()){
				if(rsidGenotype.get(key).equals("?"))
					count_missing++;
			}
			return count_missing;
		}
		
		
		
		
		
		public static void updateUserData(String caseControl, String userFilePath){
			
			// holds the rsid and the genotype value for a file - It is later used to write the file's details as an instance of the datafile
			HashMap<String, String> rsidGenotype = new HashMap<String, String>();
				
			// start time for the entire process
			long totstarttime= System.nanoTime();
			
			// String builder that holds the userID values to be updated into allSNPsFinal.arff file
			StringBuilder userIDValues = new StringBuilder();
			
			// Walk through the user files in the folder
			try {
				Files.walk(Paths.get(userFilePath)).forEach(filePath -> {
					
					long starttime = System.nanoTime();	
					updateUserMap(rsidGenotype);	
					if (Files.isRegularFile(filePath)) {
						System.out.println(filePath);	
						
						BufferedReader br1;
						
						try {
							br1 = Files.newBufferedReader(filePath);
							String strLine;
							int lineCount=0;
									
							// Read user File line by line and update hashmap with rsid and genotype data 
							while((strLine = br1.readLine()) != null){
								if(lineCount <15 )
									lineCount++;
								else{
									String[] splitLine = strLine.split("\\s+");
									if(rsidGenotype.containsKey(splitLine[0])){
										StringBuilder splLine3Val = new StringBuilder(splitLine[3]);
										if(splitLine[3].compareToIgnoreCase(splLine3Val.reverse().toString()) < 0)
											splLine3Val = splLine3Val.reverse();
										if(snpAlleleCount.get(splitLine[0]).containsKey(splLine3Val.toString())){
											rsidGenotype.put(splitLine[0], splLine3Val.toString());
										}
									}
								}	
							}
							br1.close();
						} 
						catch (Exception e) {
							// TODO Auto-generated catch block
							System.out.println("Error in reading user files");
							e.printStackTrace();
						}
						
					}
					int miss_thresh = (int)(THRESHOLD_MISSING * snpAlleleCount.size());
					System.out.println("Threshold Missing is "+ miss_thresh);
					if(countMissingValues(rsidGenotype) < miss_thresh){
					// Read each attribute name from allSNPsFinal.arff file and update the user data at the corresponding location
					FileInputStream fstream2;
					try {
						fstream2 = new FileInputStream(PRINT_FINAL);
					
						// Get the object of DataInput Stream
						DataInputStream dis2 = new DataInputStream(fstream2);
						BufferedReader br2 = new BufferedReader(new InputStreamReader(dis2));
					
						// prepare to write to a file 
						FileWriter fwriter = new FileWriter(PRINT_FINAL, true);
						BufferedWriter bw = new BufferedWriter(fwriter);
					
						// Get the userID from the filename
						String pName = filePath.getFileName().toString();
						int ind = pName.indexOf('_')+1;
						System.out.println(ind);
						StringBuilder geneSeq = new StringBuilder("'").append(pName.substring(0, (pName.indexOf('_', ind)))).append("'"); 
						//userIDValues.append("'").append(geneSeq.toString()).append("',");
						userIDValues.append(geneSeq).append(",");			
												
						boolean cont = true;
					
						// Read File line by line for rsid and populate the SNP values for a user from the hashmap
						String readRSID;
						br2.readLine(); // 1st line contains @RELATION
						br2.readLine();// 2nd Line is contains the userID attribute
						
						while(((readRSID = br2.readLine()) != null) && cont) {
							String[] splitLine = readRSID.split("\\s+");
						 
							if(splitLine[0].equals("@attribute") ){
							   									
									geneSeq.append(",").append(rsidGenotype.get(splitLine[1]));
								}
							else{
								System.out.println("the line is"+ readRSID);
								cont = false;
							}
						}
					
						// Write the user specific sequence 
						bw.newLine();
						geneSeq.append(",").append(caseControl);
						bw.write(geneSeq.toString());
				    
						br2.close();
						bw.flush();
						bw.close();
					}
					catch (Exception e) {
						// TODO Auto-generated catch block
						System.out.println("error while handling pruned document");
						e.printStackTrace();
					}
				}
					
					long elapsedtime = System.nanoTime()-starttime;
					rsidGenotype.clear();
					System.out.println(elapsedtime);
					});
					
					// Populate the user_ID attribute value
				 if(caseControl.equals("Case")){
					 System.out.println("casecontrol is case");
					 userIDValues.append(" userIdLine");
				 }else{
					 userIDValues.append("}"); 
				 }
					String pathString = "C:/Users/skatla/workspace/prepDataFile/"+PRINT_FINAL;
					Path currFile = Paths.get(pathString);
					try (Stream<String> lines = Files.lines(currFile)) {
						   List<String> replaced = lines.map(userList-> userList.replaceFirst("userIdLine", userIDValues.toString()))
						       .collect(Collectors.toList());
						   Files.write(currFile, replaced);
						}
								
			} catch (Exception e) {
				// TODO Auto-generated catch block
				System.out.println(" error while walking through the files");
				e.printStackTrace();
				}
			long totelapsedtime = System.nanoTime()-totstarttime;
			System.out.println(totelapsedtime);
		 }
		
		
		
		
		
		
		public static void main(String[] args) {
			// TODO Auto-generated method stub
			// Walk through the files to update the HashMaps
				
			try{
				Files.walk(Paths.get(READ_FROM_PATH)).forEach(filePath -> {
				if (Files.isRegularFile(filePath)) {
					BufferedReader br;
					try {
						br = Files.newBufferedReader(filePath);
						System.out.println(filePath);
						String strLine;
						int lineCount = 0;
						
						while((strLine = br.readLine()) != null){
							
							if(lineCount <15)
								lineCount++;
							else{
								String[] splitLine = strLine.split("\\s+");
								StringBuilder splLine3 = new StringBuilder();
								
								if(splitLine[3].equals("--") || splitLine[3].equals("0") || (splitLine[3].length() != 2) || (splitLine[3].equals("\\s+"))){
									splLine3.append("?");
								}							
								else{ 
									/* ensure that XY and YX are encoded as XY */
									splLine3.append(splitLine[3]);
									if(splitLine[3].compareToIgnoreCase(splLine3.reverse().toString()) < 0)
										splLine3 = splLine3.reverse();
								}
								
								TreeMap<String, Integer> innermap;
								if( snpAlleleCount.containsKey(splitLine[0])){
									Integer val;
									innermap = new TreeMap<String, Integer>();
									innermap = snpAlleleCount.get(splitLine[0]);
									
									if((val = innermap.get(splLine3.toString())) != null){
										val = val+1;
										innermap.put(splLine3.toString(), val);
									}
									else{
										innermap.put(splLine3.toString(), 1);
									}
									
									snpAlleleCount.put(splitLine[0], innermap);
								}
								else{
									innermap = new TreeMap<String, Integer>();
									innermap.put(splLine3.toString(), 1);
									snpAlleleCount.put(splitLine[0], innermap);
								}
							}
						}
						br.close();
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						e.getMessage();
						System.out.println("Error while reading the file");
					}
				}
				});
			}catch(Exception e){
				e.printStackTrace();
				e.getMessage();
				System.out.println( "Error while walking the files");
			}
					
			printsnpAlleleCount("allSNPs_June14.txt");
			
			// Select the valid SNPs - Remove outliers, discard SNPs with single alleles and populate all the alleles for double allele SNPs
			preProcessSNPAlleleCount();
			printsnpAlleleCount("preProcess_Completed.txt");
			
			// Adjust the snpValues
			Iterator<Map.Entry<String, TreeMap<String, Integer>>> it = snpAlleleCount.entrySet().iterator();
			while(it.hasNext()){
				Map.Entry<String, TreeMap<String, Integer>> ent = it.next();
				
				if(ent.getValue().size() == 2){
					TreeMap<String, Integer> innerMap1 = new TreeMap<String, Integer>();
					innerMap1 = ent.getValue();
					adjustSNPValues(innerMap1, ent.getKey());
				}
			}
					
			printsnpAlleleCount("snpValuesAdjusted.txt");
					
			// Calculate MAF
			Iterator<Map.Entry<String, TreeMap<String, Integer>>> it3 = snpAlleleCount.entrySet().iterator();
			while(it3.hasNext()){
				Map.Entry<String, TreeMap<String, Integer>> ent = it3.next();
				calcMAF(ent.getValue(), ent.getKey());
			}
			
			printsnpAlleleCount("MAF_Calculated.txt");
					
			// Remove rsid if MAF < THRESHOLD_MAF
			Iterator<Map.Entry<String, TreeMap<String, Integer>>> it4 = snpAlleleCount.entrySet().iterator();
			while(it4.hasNext()){
				Map.Entry<String, TreeMap<String, Integer>> ent = it4.next();
				if(SNPCount.get(ent.getKey()) < THRESHOLD_MAF)
					it4.remove();
			}
					
			// Write all the values to the file
			printsnpAlleleCount("final_MAFadjusted.txt");
			
			printSubjectsAndSNPs(PRINT_FINAL);
			
		}
}
