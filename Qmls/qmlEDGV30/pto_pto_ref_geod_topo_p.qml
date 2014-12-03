<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="referencialaltim">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Imbituba" value="1"/>
        <value key="Torres" value="3"/>
        <value key="Santana" value="5"/>
        <value key="Outra refer�ncia" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tiporef">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Planim�trico" value="1"/>
        <value key="Planialtim�trico" value="22"/>
        <value key="Altim�trico" value="23"/>
        <value key="Gravim�trico" value="24"/>
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
    <edittype widgetv2type="ValueMap" name="proximidade">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Isolada" value="3"/>
        <value key="Adjacente" value="14"/>
        <value key="Coincidente" value="15"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoptorefgeodtopo">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Ponto astron�mico - PA" value="1"/>
        <value key="Esta��o gravim�trica - EG" value="2"/>
        <value key="Refer�ncia de n�vel - RN" value="3"/>
        <value key="V�rtice de triangula��o - VT" value="4"/>
        <value key="Ponto barom�trico - B" value="5"/>
        <value key="Esta��o de poligonal - EP" value="6"/>
        <value key="Ponto de sat�lite - SAT" value="7"/>
        <value key="Ponto trigonom�trico - RV" value="8"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="redereferencia">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Nacional" value="1"/>
        <value key="Privada" value="2"/>
        <value key="Estadual" value="14"/>
        <value key="Municipal" value="15"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="referencialgrav">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Potsdam1930" value="2"/>
        <value key="IGSN71" value="3"/>
        <value key="Absoluto" value="4"/>
        <value key="Local" value="5"/>
        <value key="Desconhecido" value="95"/>
        <value key="N�o aplic�vel" value="97"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="situacaomarco">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o constru�do" value="1"/>
        <value key="N�o visitado" value="2"/>
        <value key="Bom" value="3"/>
        <value key="Destru�do" value="4"/>
        <value key="Destru�do sem chapa" value="5"/>
        <value key="Destru�do com chapa danificada" value="6"/>
        <value key="N�o encontrado" value="7"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>