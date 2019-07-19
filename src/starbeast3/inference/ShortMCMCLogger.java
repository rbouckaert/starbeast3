package starbeast3.inference;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

import beast.core.Logger;
import beast.core.util.Log;

public class ShortMCMCLogger extends Logger {

	@Override
	protected boolean openLogFile() throws IOException {
		String fileName = fileNameInput.get();
		final File file = new File(fileName);
		if (file.exists()) {
			if (mode == LOGMODE.compound) {
				// first find the sample nr offset
				final BufferedReader fin = new BufferedReader(new FileReader(fileName));
				String str = null;
				while (fin.ready()) {
					str = fin.readLine();
				}
				fin.close();
				assert str != null;
				final long sampleOffset = Long.parseLong(str.split("\\s")[0]);
//				if (Logger.sampleOffset > 0 && sampleOffset != Logger.sampleOffset) {
//					throw new RuntimeException("Error 400: Cannot resume: log files do not end in same sample number");
//				}
//				Logger.sampleOffset = sampleOffset;
				// open the file for appending
				final FileOutputStream out2 = new FileOutputStream(fileName, true);
				m_out = new PrintStream(out2);
			} else {
				// it is a tree logger, we may need to get rid of the last line!

				// back up file in case something goes wrong (e.g. an out of
				// memory error occurs)
				final File treeFileBackup = new File(fileName);

				// final boolean ok = treeFileBackup.renameTo(new File(fileName
				// + ".bu")); assert ok;
				Files.move(treeFileBackup.toPath(), new File(fileName + ".bu").toPath(),
						StandardCopyOption.ATOMIC_MOVE);
				// open the file and write back all but the last line
				final BufferedReader fin = new BufferedReader(new FileReader(fileName + ".bu"));

				final FileOutputStream out2 = new FileOutputStream(fileName);
				m_out = new PrintStream(out2);

				// final StringBuilder buf = new StringBuilder();
				String strLast = null;
				// String str = fin.readLine();
				boolean endSeen = false;
				while (fin.ready()) {
					if (endSeen) {
						m_out.println("End;");
						endSeen = false;
					}
					final String str = fin.readLine();
					if (!str.equals("End;")) {
						m_out.println(str);
						strLast = str;
					} else {
						endSeen = true;
					}
				}
				fin.close();

				// determine number of the last sample
				if (strLast == null) {
					// empty log file?
					throw new RuntimeException("Error 402: empty tree log file " + fileName
							+ "? (check if there is a back up file " + fileName + ".bu)");
				}
				final String str = strLast.split("\\s+")[1];
				final long sampleOffset = Long.parseLong(str.substring(6));
//				if (Logger.sampleOffset > 0 && sampleOffset != Logger.sampleOffset) {
//					// final boolean ok1 = treeFileBackup.renameTo(new
//					// File(fileName)); assert ok1;
//					Files.move(treeFileBackup.toPath(), new File(fileName).toPath(), StandardCopyOption.ATOMIC_MOVE);
//					throw new RuntimeException("Error 401: Cannot resume: log files do not end in same sample number");
//				}
//				Logger.sampleOffset = sampleOffset;
				// it is safe to remove the backup file now
				new File(fileName + ".bu").delete();
			}
			// Log.info.println("Appending file " + fileName);
			return false;
		} else {
			m_out = new PrintStream(fileName);
			//Log.warning.println("WARNING: Resuming, but file " + fileName
			//		+ " does not exist yet (perhaps the seed number is not the same as before?).");
			Log.info.println("Writing new file " + fileName);
			return true;
		}
	} // openLogFile
}
