package au.edu.uq.imb.memesuite.db;

import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * A cismapper database entry.
 */
public class CismapperDB implements JsonWr.JsonValue {
  private long listingId;
  private String name;
  private String description;
  private String genome_release;
  private String rna_source;
  private String tissues;
  private String histone_root;
  private String histone_names;
  private String max_link_distances;
  private String expression_root;
  private String expression_file_type;
  private String annotation_file_name;
  private String annotation_type;
  private String transcript_types;

  public CismapperDB(
    long listingId,
    String name,
    String description,
    String genome_release,
    String rna_source,
    String tissues,
    String histone_root,
    String histone_names,
    String max_link_distances,
    String expression_root,
    String expression_file_type,
    String annotation_file_name,
    String annotation_type,
    String transcript_types
  ) {
    this.listingId = listingId;
    this.name = name;
    this.description = description;
    this.genome_release = genome_release;
    this.rna_source = rna_source;
    this.tissues = tissues;
    this.histone_root = histone_root;
    this.histone_names = histone_names;
    this.max_link_distances = max_link_distances;
    this.expression_root = expression_root;
    this.expression_file_type = expression_file_type;
    this.annotation_file_name = annotation_file_name;
    this.annotation_type = annotation_type;
    this.transcript_types = transcript_types;
  }

  public long getListingId() {
    return listingId;
  }
  public String getName() {
    return name;
  }
  public String getDescription() {
    return description;
  }
  public String getGenomeRelease() {
    return genome_release;
  }
  public String getRnaSource() {
    return rna_source;
  }
  public String getTissues() {
    return tissues;
  }
  public String getHistoneRoot() {
    return histone_root;
  }
  public String getHistoneNames() {
    return histone_names;
  }
  public String getMaxLinkDistances() {
    return max_link_distances;
  }
  public String getExpressionRoot() {
    return expression_root;
  }
  public String getExpressionFileType() {
    return expression_file_type;
  }
  public String getAnnotationFileName() {
    return annotation_file_name;
  }
  public String getAnnotationType() {
    return annotation_type;
  }
  public String getTranscriptTypes() {
    return transcript_types;
  }

  @Override
  public void outputJson(JsonWr out) throws IOException {
    out.startObject();
    out.property("name", getName());
    out.property("description", getDescription());
    out.property("genome_release", getGenomeRelease());
    out.property("tissues", getTissues());
    out.property("histone_names", getHistoneNames());
    out.property("max_link_distances", getMaxLinkDistances());
    out.property("expression_file_type", getExpressionFileType());
    out.endObject();
  }
}
