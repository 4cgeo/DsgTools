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
        <value key="Não" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipooutlimofic">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Zona contígua" value="1"/>
        <value key="Zona econômica exclusiva" value="2"/>
        <value key="Plataforma continental jurídica" value="3"/>
        <value key="Mar territorial" value="4"/>
        <value key="Lateral marítima" value="5"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="referenciallegal">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Arruamento" value="1"/>
        <value key="Costa visível da carta" value="2"/>
        <value key="Cumeada" value="3"/>
        <value key="Limite de massa dŽágua" value="4"/>
        <value key="Linha seca" value="5"/>
        <value key="Massa dŽágua" value="6"/>
        <value key="Trecho de drenagem" value="7"/>
        <value key="Trecho ferroviário" value="8"/>
        <value key="Trecho rodoviário" value="9"/>
        <value key="Não identificado" value="10"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>