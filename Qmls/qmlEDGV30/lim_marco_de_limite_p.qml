<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipomarcolim">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Municipal" value="3"/>
        <value key="Estadual" value="23"/>
        <value key="Internacional secund�rio" value="24"/>
        <value key="Internacional de refer�ncia" value="25"/>
        <value key="Internacional principal" value="26"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="sistemageodesico">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Astro Chu�" value="1"/>
        <value key="C�rrego Alegre" value="2"/>
        <value key="SAD-69" value="3"/>
        <value key="WGS-84" value="5"/>
        <value key="SIRGAS2000" value="6"/>
        <value key="Outra refer�ncia" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="referencialaltim">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Imbituba" value="1"/>
        <value key="Torres" value="3"/>
        <value key="Santana" value="5"/>
        <value key="Outra refer�ncia" value="99"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>