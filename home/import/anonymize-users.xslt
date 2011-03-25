<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
       xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
       version='1.0'
       >
<xsl:output method="xml" indent="yes"/>

<xsl:template match="/">
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="*|@*">
<xsl:copy>
<xsl:apply-templates select="*|@*|text()"/>
</xsl:copy>
</xsl:template>

<xsl:template match="OpenId">
<OpenId>
       <xsl:text>http://nowhere.com</xsl:text>
</OpenId>
</xsl:template>

<xsl:template match="Birthday">
<Birthday>
       <xsl:text>1900-01-01</xsl:text>
</Birthday>
</xsl:template>

<xsl:template match="Email">
<Email>
       <xsl:text>no.me@nowhere.com</xsl:text>
</Email>
</xsl:template>


</xsl:stylesheet>