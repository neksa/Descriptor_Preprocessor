package au.edu.uq.imb.memesuite.servlet.util;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import java.sql.SQLException;
import java.util.EnumSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.paramRequire;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CISMAPPER_DB_KEY;

/**
 * Data entry component for Cismapper databases
 */
public class ComponentCismapper extends PageComponent {
  private ServletContext context;
  private HTMLTemplate tmplCismapper;
  private HTMLTemplate tmplCategory;
  private HTMLTemplate tmplListing;
  private String prefix;
  private HTMLTemplate title;
  private HTMLTemplate subtitle;


  private static final Pattern DB_ID_PATTERN = Pattern.compile("^\\d+$");
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.component.cismapper");

  public ComponentCismapper(ServletContext context, HTMLTemplate info) throws ServletException {
    this.context = context;
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplCismapper = cache.loadAndCache("/WEB-INF/templates/component_cismapper.tmpl");
    prefix = getText(info, "prefix", "motifs");
    title = getTemplate(info, "title", null);
    subtitle = getTemplate(info, "subtitle", null);
    tmplCategory = tmplCismapper.getSubtemplate("component").getSubtemplate("category");
    tmplListing = tmplCategory.getSubtemplate("listing");
  }

  @Override
  public HTMLSub getComponent() {
    HTMLSub cismapper = tmplCismapper.getSubtemplate("component").toSub();
    cismapper.set("prefix", prefix);
    if (title != null) cismapper.set("title", title);
    if (subtitle != null) cismapper.set("subtitle", subtitle);
    CismapperDBList db = (CismapperDBList)context.getAttribute(CISMAPPER_DB_KEY);
    if (db == null) {
      logger.log(Level.SEVERE, "No Cismapper database is intalled");
      cismapper.empty("category");
      return cismapper;
    }
    try {
      cismapper.set("category", new AllCategories(db));
    } catch (SQLException e) {
      logger.log(Level.SEVERE, "Error querying Cismapper categories", e);
      cismapper.empty("category");
    }
    return cismapper;
  }

  @Override
  public HTMLSub getHelp() {
    return this.tmplCismapper.getSubtemplate("help").toSub();
  }

  public CismapperDB getCismapperPanel(HttpServletRequest request) throws ServletException {
    // determine the source
    String source = paramRequire(request, prefix + "_source");
    Matcher m = DB_ID_PATTERN.matcher(source);
    if (!m.matches()) {
      throw new ServletException("Parameter " + prefix + "_source had a " +
          "value that did not match any of the allowed values.");
    }
    long dbId;
    try {
      dbId = Long.parseLong(source, 10);
    } catch (NumberFormatException e) {
      throw new ServletException(e);
    }
    CismapperDBList db = (CismapperDBList)context.getAttribute(CISMAPPER_DB_KEY);
    if (db == null) {
      throw new ServletException("Unable to access the Cismapper database.");
    }
    try {
      return db.getCismapperListing(dbId);
    } catch (SQLException e) {
      throw new ServletException(e);
    }
  } // getCismapperPanel

  private class AllCategories extends HTMLSubGenerator<Category> {
    private CismapperDBList db;
    private AllCategories(CismapperDBList db) throws SQLException {
      super(db.getCategories(false, EnumSet.allOf(AlphStd.class)));
      this.db = db;
    }
    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = tmplCategory.toSub();
      out.set("name", item.getName());
      try {
        out.set("listing", new AllListingsOfCategory(db, item.getID()));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error querying Cismapper listings", e);
        out.empty("listing");
      }
      return out;
    }
  }

  private class AllListingsOfCategory extends HTMLSubGenerator<Listing> {
    private AllListingsOfCategory(CismapperDBList db, long id) throws SQLException {
      super(db.getListings(id));
    }

    @Override
    protected HTMLSub transform(Listing item) {
      HTMLSub out = tmplListing.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      return out;
    }
  }
} //class ComponentCismapper
