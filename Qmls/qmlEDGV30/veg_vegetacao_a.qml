<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="TextEdit" name="id_area_verde"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoveg">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Vegeta��o cultivada" value="1"/>
        <value key="Floresta" value="2"/>
        <value key="Vegeta��o de mangue" value="3"/>
        <value key="Ref�gio ecol�gico" value="4"/>
        <value key="Campinarana" value="5"/>
        <value key="Cerrado" value="6"/>
        <value key="Vegeta��o de restinga" value="7"/>
        <value key="Estepe" value="8"/>
        <value key="Vegeta��o de brejo ou p�ntano" value="9"/>
        <value key="Caatinga" value="10"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="classificacaoporte">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Rasteira" value="2"/>
        <value key="Herb�cea" value="3"/>
        <value key="Arb�rea" value="4"/>
        <value key="Arbustiva" value="5"/>
        <value key="Desconhecida" value="95"/>
        <value key="Mista" value="97"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>