package cophesim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;

public class ConvertCophesimToPLINK {
	static Path phenoFile, pedFile;
	static String outPrefix = "PLINK";
	static double minCasePercent = 0.3;
	static double maxCasePercent = 0.7;
	static int nSNPs, nSamples, nCausal;
	static int repetitions = 10;
	static Path outDir = Paths.get("CophesimPLINK");
	static Path plinkParamFile, effectsFile;
	public static void main(String[] args) {
		if(args.length < 1){
			printHelp();
			return;
		}
		parseArgs(args);
		if(!Files.exists(outDir)){
			try {
				Files.createDirectories(outDir);
			} catch (IOException e) {
				System.out.println("Error while creating outDir");
				System.out.println(e.getMessage());
				return;
			}
		}
		effectsFile = Paths.get(outDir + "/effects.txt");
		plinkParamFile = Paths.get(outDir + "/plinkParam");
		try{
			Files.write(plinkParamFile, (nSNPs + " null 0.00 1.00 1.00 1.00").getBytes());
		}
		catch(IOException e){
			System.out.println("Error while writing plinkParamFile");
			System.out.println(e.getMessage());
			return;
		}
		int trialNumber = 0;
		double casePercent;
		do{		
			int innerTrial = 0;
			try{
				execPLINKSim();
			}	
			catch(IOException e){
				System.out.println("Error while simulating data");
				System.out.println(e.getMessage());
				return;
			}
			do{
			try{
				List<Integer> pos = ThreadLocalRandom.current().ints(0, nSNPs).limit(nCausal).boxed().collect(Collectors.toList());
				List<Integer> pos2 = ThreadLocalRandom.current().ints(0, nSNPs).limit(nCausal).boxed().collect(Collectors.toList());
				List<Double> effect = ThreadLocalRandom.current().doubles(-2, 2).limit(nCausal).boxed().sorted().collect(Collectors.toList());
				Files.deleteIfExists(effectsFile);
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(effectsFile.toFile())));
				for(int i = 0; i < nCausal; i++){
					//pw.println(pos.get(i) + ":" + effect.get(i));
					pw.println(pos.get(i) + "," + pos2.get(i) + "," + effect.get(i));
				}
				pw.close();
				casePercent = execCophesim();
				innerTrial++;
				System.out.println("\tInnerTrial " + innerTrial + " = " + casePercent);
			}
			catch(IOException e){
				System.out.println("Error while simulating data");
				System.out.println(e.getMessage());
				return;
			}
			}
			while(innerTrial < repetitions && (casePercent < minCasePercent || casePercent > maxCasePercent));
			System.out.println("Trial " + trialNumber + " = " + casePercent);
			trialNumber++;
		}
		while(casePercent < minCasePercent || casePercent > maxCasePercent);
		phenoFile = Paths.get(outDir.toAbsolutePath().toString() + "/cophesimOut_pheno_bin.txt.pheno");
		pedFile = Paths.get(outDir.toAbsolutePath().toString() + "/cophesimOut_pheno_bin.txt.ped");
		try{
			Files.copy(Paths.get(outDir.toAbsolutePath().toString() + "/cophesimOut_pheno_bin.txt.map"),Paths.get(outDir.toAbsolutePath().toString() + "/PLINK.map"), StandardCopyOption.REPLACE_EXISTING);
		}
		catch(IOException e){
			System.out.println("Error while copying mapFile");
			System.out.println(e.getMessage());
			return;
		}
		Map<Integer,Integer> phenotypeMap;
		try{
			phenotypeMap = readPhenotypes();
		}
		catch(IOException e){
			System.out.println("Error while reading phenoFile");
			System.out.println(e.getMessage());
			return;
		}	
		try{
			printPedFile(phenotypeMap);
		}
		catch(IOException e){
			System.out.println("Error while writing to pedFile");
			System.out.println(e.getMessage());
			return;
		}
	}
	private static void execPLINKSim() throws IOException{
		List<String> commandList = new ArrayList<>();
		commandList.add("plink");
		commandList.add("--simulate-ncases");
		commandList.add(Integer.toString(nSamples / 2));
		commandList.add("--simulate-ncontrols");
		commandList.add(Integer.toString(nSamples / 2));
		commandList.add("--simulate");
		commandList.add(plinkParamFile.toAbsolutePath().toString());
		commandList.add("--out");
		commandList.add("plinkSimu");
		commandList.add("--make-bed");
		ProcessBuilder pbPlink = new ProcessBuilder(commandList);
		pbPlink.directory(outDir.toFile());
		try {
			execProc(pbPlink, true);
		} catch (IOException | InterruptedException e) {
			throw new IOException("Error while using PLINK\n" + e.getMessage());
		}
	}
	private static double execCophesim() throws IOException{
		List<String> commandList = new ArrayList<>();
		commandList.clear();
		commandList.add("python");
		commandList.add("/Home/felix.heinrich/git/AZIFI/src/cophesim/Scripts/cophesim.py");
		commandList.add("-i");
		commandList.add("plinkSimu");
		commandList.add("-o");
		commandList.add("cophesimOut");
		commandList.add("-epi");
		commandList.add(effectsFile.toAbsolutePath().toString());
		ProcessBuilder pbCophesim = new ProcessBuilder(commandList);
		pbCophesim.directory(outDir.toFile());
		try {
			execProc(pbCophesim, true);
		} catch (IOException | InterruptedException e) {
			throw new IOException("Error while using Cophesim\n" + e.getMessage());
		}
		long nCase = Files.readAllLines(Paths.get(outDir.toAbsolutePath().toString() + "/cophesimOut_pheno_bin.txt.pheno")).stream()
			.map(line -> line.split(" ")[2]).filter(number -> number.equals("2")).collect(Collectors.counting());
		return (double) nCase / nSamples;
	}
	private static Map<Integer, Integer> readPhenotypes() throws IOException{
		Map<Integer,Integer> phenotypeMap = new HashMap<>();
		List<String> lines = Files.readAllLines(phenoFile);
		for(int i = 0; i < lines.size(); i++){
			String[] tmp = lines.get(i).split("\\s+");
			phenotypeMap.put(Integer.parseInt(tmp[0]), Integer.parseInt(tmp[2]));
		}
		return phenotypeMap;
	}
	private static void printPedFile(Map<Integer,Integer> phenotypeMap) throws IOException{
		Path outFile = Paths.get(outDir.toAbsolutePath().toString() + "/PLINK.ped");
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outFile.toFile())));
		List<String> lines = Files.readAllLines(pedFile);
		for(int i = 0; i < lines.size(); i++){
			String[] tmp = lines.get(i).split("\\s+");
			int id = Integer.parseInt(tmp[0]);
			pw.print("Test" + id + " Test" + tmp[1] + " " + tmp[2] + " " + tmp[3] + " " + tmp[4] + " " + phenotypeMap.get(id));
			for(int j = 6; j < tmp.length; j++){
				pw.print(" " + tmp[j]);
			}
			pw.println();
		}
		pw.close();
	}
	private static void execProc(ProcessBuilder pb, boolean redirectError) throws IOException, InterruptedException{
		if(redirectError){
			pb.redirectErrorStream(true);
		}
		else{
			pb.redirectErrorStream(false);
		}
		Process proc = pb.start();
		String line;
		BufferedReader procIn = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		while((line = procIn.readLine()) != null) {
			System.out.println(line);
		}
		procIn.close();
		proc.waitFor();
	}
	private static void parseArgs(String[] args){
		for(int i = 0; i < args.length - 3; i++){
			switch(args[i]){
			case "-out" : 
				outDir = Paths.get(args[i+1]);
				i++;
				break;
			}
		}
		nSNPs = Integer.parseInt(args[args.length-3]);
		nSamples = Integer.parseInt(args[args.length-2]);
		nCausal = Integer.parseInt(args[args.length-1]);
	}
	private static void printHelp(){
		System.out.println("Usage: java -jar ConvertCophesimToPLINK.jar {Options} nSNPs nSamples nCausal");
		System.out.println("Options:");
		System.out.println("-out dir\t name of output directory (default CophesimPLINK)");
	}
}
