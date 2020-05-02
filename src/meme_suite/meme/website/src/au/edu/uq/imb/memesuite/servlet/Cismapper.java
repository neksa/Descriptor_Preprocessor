package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.LociDataSource;
import au.edu.uq.imb.memesuite.db.CismapperDB;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.activation.DataSource;
import javax.servlet.*;
import javax.servlet.http.*;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;


public class Cismapper extends SubmitJob<Cismapper.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentLoci loci;
  private ComponentCismapper cismapperPanel;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advBtn;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.cismapper");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public LociDataSource loci;
    public CismapperDB cismapperPanel;
    String rna_source;
    String tissues;
    String histone_root;
    String histone_names;
    String max_link_distances;
    String expression_root;
    String expression_file_type;
    String annotation_file_name ;
    String annotation_type;
    String transcript_types;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      out.property("loci", loci);
      out.property("cismapperPanel", cismapperPanel);
      out.property("rna_source", rna_source);
      out.property("tissues", tissues);
      out.property("histone_root", histone_root);
      out.property("histone_names", histone_names);
      out.property("max_link_distances", max_link_distances);
      out.property("expression_root", expression_root);
      out.property("expression_file_type", expression_file_type);
      out.property("annotation_file_name", annotation_file_name);
      out.property("annotation_type", annotation_type);
      out.property("transcript_types", transcript_types);
      out.endObject();
    }

    @Override
    public String email() {
      return email;
    }
  
    @Override
    public String description() {
      return description;
    }

    @Override
    public String emailTemplate() {
      return tmplVerify.getSubtemplate("message").toString();
    }

    @Override
    public String cmd() {
  //    cismapper_webservice [options] <locus_file> <rna_source> 
  //
  //      Options:
  //
  //	-tissues <tissues>		comma-separated list (no spaces) of tissue names 
  //	-histone-root <hrd>		histone root directory 
  //	-histone-names <hnames>		comma-separated list (no spaces) of histone names
  //	-max-link-distances <mlds>	comma-separated list of maximum distances between an RE and its target
  //	-expression-root <erd>		expression root directory
  //	-expression-file-type <eft>	file type of expression files
  //	-annotation-file-name <afile>	annotation file name
  //	-annotation-type <atype>	type of annotation [Gencode|RefSeq]
  //	-transcript-types <ttypes>	types of transcript to use from annotation file
  //	-min-feature-count <mfc>	only consider links where there is both histone and
  //					expression data for at least this many tissues: default: 7
  //	-min-max-expression <mme>	maximum expression of a target must be at least <mme> for
  //					the target to be included in the map; default: 2
  //	-max-html-score <mhs>		only include links with this score or better in the HTML
  //					default: 0.05
  //
  //      Files present in the server cismapper databases can be specified by appending 'db/'
  //      to the file name.
      StringBuilder args = new StringBuilder();
      addArgs(args, loci.getName());
      addArgs(args, cismapperPanel.getRnaSource());
      addArgs(args, "-tissues", cismapperPanel.getTissues());
      addArgs(args, "-histone-root", "db/"+cismapperPanel.getHistoneRoot());
      addArgs(args, "-histone-names", cismapperPanel.getHistoneNames());
      addArgs(args, "-max-link-distances", cismapperPanel.getMaxLinkDistances());
      addArgs(args, "-expression-root", "db/"+cismapperPanel.getExpressionRoot());
      addArgs(args, "-expression-file-type", cismapperPanel.getExpressionFileType());
      addArgs(args, "-annotation-file-name", "db/"+cismapperPanel.getAnnotationFileName());
      addArgs(args, "-annotation-type", cismapperPanel.getAnnotationType());
      addArgs(args, "-transcript-types", cismapperPanel.getTranscriptTypes());
      addArgs(args, "-min-feature-count", "3");
      addArgs(args, "-min-max-expression", "2");
      addArgs(args, "-max-html-score", "0.05");
      addArgs(args, "-nostatus");
      addArgs(args, "-noecho");
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (loci != null) list.add(loci);
      return list;
    }
  
    @Override
    public void cleanUp() {
      if (loci != null) {
        if (!loci.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              loci.getFile());
        }
      }
    }
  }

  public Cismapper() {
    super("CISMAPPER", "CisMapper");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplMain = cache.loadAndCache("/WEB-INF/templates/cismapper.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/cismapper_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    loci = new ComponentLoci(context, tmplMain.getSubtemplate("loci"));
    cismapperPanel = new ComponentCismapper(context, tmplMain.getSubtemplate("cismapper_panel"));
    jobDetails = new ComponentJobDetails(cache);
    advBtn = new ComponentAdvancedOptions(cache);
    submitReset = new ComponentSubmitReset(cache, jobTable.getCount(), jobTable.getDuration());
    footer = new ComponentFooter(cache, msp);
  }

  @Override
  public String title() {
    return tmplVerify.getSubtemplate("title").toString();
  }

  @Override
  public String subtitle() {
    return tmplVerify.getSubtemplate("subtitle").toString();
  }

  @Override
  public String logoPath() {
    return tmplVerify.getSubtemplate("logo").toString();
  }

  @Override
  public String logoAltText() {
    return tmplVerify.getSubtemplate("alt").toString();
  }

  @Override
  protected void displayForm(HttpServletRequest request, HttpServletResponse response, long quotaMinWait) throws IOException {
    HTMLSub main = tmplMain.toSub();
    main.set("help", new HTMLSub[]{header.getHelp(), loci.getHelp(),
        cismapperPanel.getHelp(), jobDetails.getHelp(), advBtn.getHelp(),
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("loci", loci.getComponent());
    main.set("cismapper_panel", cismapperPanel.getComponent());
    main.set("job_details", jobDetails.getComponent());
    // main.set("advanced_options", advBtn.getComponent());
    main.set("submit_reset", submitReset.getComponent(quotaMinWait));
    main.set("footer", footer.getComponent());
    response.setContentType("text/html; charset=UTF-8");
    main.output(response.getWriter());
  }

  @Override
  protected Data checkParameters(FeedbackHandler feedback,
      HttpServletRequest request) throws IOException, ServletException {
    // setup default file names
    FileCoord namer = new FileCoord();
    FileCoord.Name lociName = namer.createName("loci.bed");
    namer.createName("description");
    namer.createName("uuid");
    Data data = new Data();
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // get the loci
    data.loci = (LociDataSource)loci.getLoci(lociName, request, feedback);
    // get the panel 
    data.cismapperPanel = cismapperPanel.getCismapperPanel(request);

    return data;
  } // checkParameters
}
