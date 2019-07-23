package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.servlet.util.WebUtils;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

// debugging FIXME
//import java.lang.Object;
//import java.io.Serializable;
//import org.apache.commons.lang.builder.ToStringBuilder;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CISMAPPER_DB_KEY;

/**
 * Show the databases available to the CISMAPPER program.
 */
public class ShowCismapperDBs extends HttpServlet {
  private ServletContext context;
  private HTMLTemplate template;
  private HTMLTemplate categoryTemplate;
  private HTMLTemplate listingTemplate;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web");

  public ShowCismapperDBs() { }

  @Override
  public void init() throws ServletException {
    context = this.getServletContext();
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    template = cache.loadAndCache("/WEB-INF/templates/show_cismapper_dbs.tmpl");
    categoryTemplate = template.getSubtemplate("content").getSubtemplate("category");
    listingTemplate = categoryTemplate.getSubtemplate("listing");
  }

  @Override
  public void doGet(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    if (request.getParameter("category") != null) {
      outputXmlListingsOfCategory(response, getId(request, "category"));
    } else {
      display(response);
    }
  }

  @Override
  public void doPost(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    display(response);
  }

  private void display(HttpServletResponse response)
      throws IOException, ServletException {
    CismapperDBList cismapperDBList = (CismapperDBList)context.getAttribute(CISMAPPER_DB_KEY);
    response.setContentType("text/html; charset=UTF-8");
    HTMLSub out = template.toSub();
    if (cismapperDBList != null) {
      try {
        out.getSub("content").set("category", new AllCategories(cismapperDBList));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading categories", e);
        out.empty("content");
      }
    } else {
      out.empty("content");
    }
    out.output(response.getWriter());
  }

  private long getId(HttpServletRequest request, String name) throws ServletException {
    String value = request.getParameter(name);
    if (value == null) {
      throw new ServletException("Parameter '" + name + "' was not set.");
    }
    long id;
    try {
      id = Long.parseLong(value, 10);
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter '" + name + "' is not a integer value.", e);
    }
    return id;
  }

  private void outputXmlListingsOfCategory(HttpServletResponse response, long categoryId)
      throws IOException, ServletException {
    CismapperDBList cismapperDBList = (CismapperDBList)context.getAttribute(CISMAPPER_DB_KEY);
    response.setContentType("application/xml; charset=UTF-8");
    PrintWriter out = null;
    try {
      Iterator<Listing> listingView = cismapperDBList.getListings(categoryId).iterator();
      out = response.getWriter();
      out.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
      out.println("<listings category=\"" + categoryId + "\">");
      while (listingView.hasNext()) {
        Listing listing = listingView.next();
        out.println("<l i=\"" + listing.getID() + "\" n=\"" + WebUtils.escapeForXML(listing.getName()) + "\"/>");
      }
      out.println("</listings>");
    } catch (SQLException e) {
      throw new ServletException("Error getting CISMAPPER listings", e);
    } finally {
      if (out != null) out.close();
    }
  }

  private class AllCategories extends HTMLSubGenerator<Category> {
    private CismapperDBList cismapperDBList;
    private AllCategories(CismapperDBList cismapperDBList) throws SQLException {
      super(cismapperDBList.getCategories(false, EnumSet.allOf(AlphStd.class)));
      this.cismapperDBList = cismapperDBList;
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = categoryTemplate.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      try {
        out.set("listing", new AllListingsOfCategory(cismapperDBList, item.getID()));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading CISMAPPER listings", e);
        out.empty("listing");
      }
      return out;
    }
  }

  private class AllListingsOfCategory extends HTMLSubGenerator<Listing> {
    private CismapperDBList cismapperDBList;
    private AllListingsOfCategory(CismapperDBList cismapperDBList, long id) throws SQLException{
      super(cismapperDBList.getListings(id));
      this.cismapperDBList = cismapperDBList;
    }

    @Override
    protected HTMLSub transform(Listing listing) {
      HTMLSub out = listingTemplate.toSub();
      out.set("name", listing.getName());
      out.set("description", listing.getDescription());
      try {
        CismapperDB cismapperDB = cismapperDBList.getCismapperListing(listing.getID());
        out.set("genome_release", cismapperDB.getGenomeRelease());
        out.set("tissues", cismapperDB.getTissues().replace(",", " "));
        out.set("histone_names", cismapperDB.getHistoneNames().replace(",", " "));
        out.set("max_link_distances", cismapperDB.getMaxLinkDistances().replace(",", " "));
        out.set("rna_source", cismapperDB.getRnaSource());
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error querying CISMAPPER listing",  e);
      }
      return out;
    }
  }

}
