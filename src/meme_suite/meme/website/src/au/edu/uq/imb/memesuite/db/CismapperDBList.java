package au.edu.uq.imb.memesuite.db;

import java.io.*;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.*;
import java.sql.*;

import au.edu.uq.imb.memesuite.data.AlphStd;
import org.sqlite.*;

public class CismapperDBList extends DBList {
  private File tsv;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.cismapperdb");

  public CismapperDBList(File tsv) throws ClassNotFoundException, IOException, SQLException {
    super(loadCismapperTSV(tsv), true);
    this.tsv = tsv;
  }

  /**
   * Get the source file that the database was created from.
   * @return the TSV file that was loaded to create the database.
   */
  public File getTSV() {
    return tsv;
  }

  @Override
  protected PreparedStatement prepareListingsQuery(Connection conn, long categoryId, boolean shortOnly, EnumSet<AlphStd> allowedAlphabets) throws SQLException {
    PreparedStatement ps = conn.prepareStatement(SQL.SELECT_CISMAPPER_LISTINGS_OF_CATEGORY);
    ps.setLong(1, categoryId);
    //ps.setInt(2, SQL.enumsToInt(allowedAlphabets));
    return ps;
  }

  public CismapperDB getCismapperListing(long listingId) throws SQLException {
    CismapperDB cismapperDB;
    Connection conn = null;
    PreparedStatement stmt = null;
    ResultSet rset = null;
    try {
      // create a database connection
      conn = ds.getConnection();
      stmt = conn.prepareStatement(SQL.SELECT_CISMAPPER_LISTING);
      stmt.setLong(1, listingId);
      rset = stmt.executeQuery();
      if (!rset.next()) throw new SQLException("No listing by that ID");
      String name = rset.getString(1);
      String description = rset.getString(2);
      String genome_release = rset.getString(3);
      String rna_source = rset.getString(4);
      String tissues = rset.getString(5);
      String histone_root = rset.getString(6);
      String histone_names = rset.getString(7);
      String max_link_distances = rset.getString(8);
      String expression_root = rset.getString(9);
      String expression_file_type = rset.getString(10);
      String annotation_file_name = rset.getString(11);
      String annotation_type = rset.getString(12);
      String transcript_types = rset.getString(13);
      rset.close(); rset = null;
      stmt.close(); stmt = null;
      conn.close(); conn = null;
      cismapperDB = new CismapperDB( listingId, name, description, genome_release, 
	rna_source, tissues, histone_root, histone_names, max_link_distances, 
	expression_root, expression_file_type, annotation_file_name, 
	annotation_type, transcript_types
      );
    } finally {
      if (rset != null) {
        try {
          rset.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (stmt != null) {
        try {
          stmt.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (conn != null) {
        try {
          conn.close();
        } catch (SQLException e) {
          logger.log(Level.SEVERE, "Cleanup (after error) failed to close the DB connection", e);
        }
      }
    }
    return cismapperDB;
  } // getCismapperListing

  private static File loadCismapperTSV(File tsv) throws ClassNotFoundException, IOException, SQLException {
    // load the JDBC needed
    Class.forName("org.sqlite.JDBC");
    // create the file to contain the database
    File db = File.createTempFile("cismapper_db_", ".sqlite");
    db.deleteOnExit();
    // configure the database
    SQLiteConfig config = new SQLiteConfig();
    config.enforceForeignKeys(true);
    SQLiteDataSource ds = new SQLiteDataSource(config);
    ds.setUrl("jdbc:sqlite:" + db);
    // open a connection
    Connection connection = null;
    try {
      connection = ds.getConnection();
      Statement statement = connection.createStatement();
      statement.executeUpdate(SQL.CREATE_TBL_CATEGORY);
      statement.executeUpdate(SQL.CREATE_TBL_LISTING);
      statement.executeUpdate(SQL.CREATE_TBL_CISMAPPER);
      connection.setAutoCommit(false);
      importCismapperTSV(tsv, connection);
      connection.commit();
    } finally {
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) { /* ignore */ }
      }
    }
    if (!db.setLastModified(tsv.lastModified())) {
      logger.log(Level.WARNING, "Failed to set last modified date on " + db);
    }
    logger.log(Level.INFO, "Loaded CISMAPPER TSV \"" + tsv + "\" into \"" + db + "\"");
    return db;
  } // loadCismapperTSV

  private static void importCismapperTSV(File tsv, Connection connection) throws IOException, SQLException {
    String line;
    Long categoryId = null;
    Pattern emptyLine = Pattern.compile("^\\s*$");
    Pattern commentLine = Pattern.compile("^#");
    BufferedReader in = null;
    try {
      // open the tsv file for reading
      in = new BufferedReader(new InputStreamReader(new FileInputStream(tsv), "UTF-8"));
      // create the prepared statements
      PreparedStatement pstmtCategory = connection.prepareStatement(
          SQL.INSERT_CATEGORY, Statement.RETURN_GENERATED_KEYS);
      PreparedStatement pstmtListing = connection.prepareStatement(
          SQL.INSERT_LISTING, Statement.RETURN_GENERATED_KEYS);
      PreparedStatement pstmtCismapperDB = connection.prepareStatement(
          SQL.INSERT_CISMAPPER_DB, Statement.NO_GENERATED_KEYS);
      // now read the tsv file
      while ((line = in.readLine()) != null) {
        // skip any empty lines or comments
        if (emptyLine.matcher(line).find()) continue;
        if (commentLine.matcher(line).find()) continue;
        line = line.trim();
        // check we have enough items on the line to do something
        String[] values = line.split("\\s*\t\\s*");	// split on tab
        if (values.length < 13) {
          logger.log(Level.WARNING, "TSV line has " + values.length +
              " values but expected 13.");
          continue;
        }
        // check that a name was supplied
        if (emptyLine.matcher(values[0]).find()) {
          logger.log(Level.WARNING, "TSV line has no entry for name column.");
          continue;
        }
        // test to see if we have a category or a selectable listing
	// listing
	if (values.length < 13) {
	  logger.log(Level.WARNING, "TSV line has " + values.length +
	      " values but expected 13 for a CISMAPPER listing.");
	  continue;
	}
	int i = 0;
	String name = values[i++];
	String description = values[i++];
	String genome_release = values[i++];
	String rna_source = values[i++];
	String tissues = values[i++];
	String histone_root = values[i++];
	String histone_names = values[i++];
	String max_link_distances = values[i++];
	String expression_root = values[i++];
	String expression_file_type = values[i++];
	String annotation_file_name = values[i++];
	String annotation_type = values[i++];
	String transcript_types = values[i++];
	// check we have a category to store the listing under
	if (categoryId == null) {
	  // create a dummy category with an empty name
	  pstmtCategory.setString(1, "");
	  pstmtCategory.executeUpdate();
	  // now get the category Id
	  ResultSet rs = pstmtCategory.getGeneratedKeys();
	  if (!rs.next()) throw new SQLException("Failed to get Category Id.\n");
	  categoryId = rs.getLong(1);
	  rs.close();
	}
	// now create the listing
	pstmtListing.setLong(1, categoryId);
	pstmtListing.setString(2, name);
	pstmtListing.setString(3, description);
	pstmtListing.executeUpdate();
	ResultSet rs = pstmtListing.getGeneratedKeys();
	if (!rs.next()) throw new SQLException("Failed to get Listing Id.\n");
	long listingId = rs.getLong(1);
	rs.close();
	// create the CISMAPPER primary listing
	i = 1; 
	pstmtCismapperDB.setLong(i++, listingId);
	pstmtCismapperDB.setString(i++, genome_release);
	pstmtCismapperDB.setString(i++, rna_source);
	pstmtCismapperDB.setString(i++, tissues);
	pstmtCismapperDB.setString(i++, histone_root);
	pstmtCismapperDB.setString(i++, histone_names);
	pstmtCismapperDB.setString(i++, max_link_distances);
	pstmtCismapperDB.setString(i++, expression_root);
	pstmtCismapperDB.setString(i++, expression_file_type);
	pstmtCismapperDB.setString(i++, annotation_file_name);
	pstmtCismapperDB.setString(i++, annotation_type);
	pstmtCismapperDB.setString(i++, transcript_types);
	pstmtCismapperDB.executeUpdate();
	logger.log(Level.FINE, "Loaded Cismapper Listing: " + name);
      }
      // close the prepared statements
      pstmtCategory.close();
      pstmtListing.close();
      pstmtCismapperDB.close();
    } finally {
      if (in != null) {
        try {
          in.close();
        } catch (IOException e) {
          // ignore
        }
      }
    } // try
  } // importCismapperTSV
} // class CismapperDBList 
